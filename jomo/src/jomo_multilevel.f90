! Homogeneous multilevel joint-model MCMC kernels for jomo.
! This module covers the computational core shared by jomo1rancon,
! jomo1rancat, jomo1ranmix and their chain variants.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_multilevel
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state
   use jomo_linalg, only : inverse_spd, is_spd
   use jomo_distributions, only : matrix_normal_sample, mvnormal_sample, invwishart_sample
   use jomo_latent, only : latent_dimension, initialize_latent_data, set_categorical_covariance
   use jomo_latent, only : decode_categories, update_observed_categories, impute_missing_rows
   use jomo_single_level, only : jomo1_result
   implicit none
   private

   type, public :: jomo1ran_result
      type(jomo1_result) :: level1
      real(dp), allocatable :: random_effects(:, :, :)
      real(dp), allocatable :: random_covariance(:, :)
      real(dp), allocatable :: random_covariance_mean(:, :)
      real(dp), allocatable :: random_effects_mean(:, :, :)
   end type jomo1ran_result

   public :: jomo1ran_mixed_mcmc
   public :: jomo1rancon_mcmc
   public :: compute_random_mean
   public :: random_design_row
   public :: pack_random_matrix
   public :: unpack_random_vector
   public :: mh_covariance_multilevel

contains

   subroutine jomo1ran_mixed_mcmc(rng, y_con, con_observed, y_cat, cat_observed, n_levels, x, z, cluster, n_iter, &
      prior_scale, prior_random_scale, result, beta_start, omega_start, random_start, random_cov_start)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling Gibbs, latent, and covariance draws.
      real(dp), intent(in) :: y_con(:, :) !! Continuous level-1 responses, shape n by n_con; masked entries are ignored.
      logical, intent(in) :: con_observed(:, :) !! Observation mask for continuous level-1 responses.
      integer, intent(in) :: y_cat(:, :) !! Categorical level-1 responses encoded 1..K, shape n by n_cat.
      logical, intent(in) :: cat_observed(:, :) !! Observation mask for categorical level-1 responses.
      integer, intent(in) :: n_levels(:) !! Number of levels K for each categorical level-1 response.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by qx.
      real(dp), intent(in) :: z(:, :) !! Fully observed random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based cluster label for each row; labels must be contiguous 1..G.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale(:, :) !! Residual inverse-Wishart scale matrix, p by p, used for all-continuous outcomes.
      real(dp), intent(in) :: prior_random_scale(:, :) !! Random-effect inverse-Wishart scale matrix, dimension qz*p by qz*p.
      type(jomo1ran_result), intent(out) :: result !! Final imputations, effects, covariance states, and posterior means.
      real(dp), intent(in), optional :: beta_start(:, :) !! Optional qx by p starting fixed-effect matrix; defaults to zero.
      real(dp), intent(in), optional :: omega_start(:, :) !! Optional p by p starting residual covariance; defaults to identity.
      real(dp), intent(in), optional :: random_start(:, :, :) !! Optional G by qz by p random-effect start; default zero.
      real(dp), intent(in), optional :: random_cov_start(:, :) !! Optional qz*p square random covariance start; default identity.
      integer :: n
      integer :: n_con
      integer :: n_cat
      integer :: p
      integer :: qx
      integer :: qz
      integer :: g_count
      integer :: d
      integer :: iter
      integer :: g
      integer :: i
      integer :: j
      integer :: info
      integer :: failures
      integer :: total_failures
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: omega(:, :)
      real(dp), allocatable :: u(:, :, :)
      real(dp), allocatable :: covu(:, :)
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: xtx_inv(:, :)
      real(dp), allocatable :: bhat(:, :)
      real(dp), allocatable :: random_mean(:, :)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: target(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: scale(:, :)
      real(dp), allocatable :: covu_precision(:, :)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: post_precision(:, :)
      real(dp), allocatable :: post_cov(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: post_mean(:)
      real(dp), allocatable :: uvec(:)
      real(dp), allocatable :: hmat(:, :)
      logical, allocatable :: observed(:, :)
      logical, allocatable :: fixed(:, :)

      n = size(x, 1)
      if (size(z, 1) /= n .or. size(cluster) /= n) error stop "jomo1ran_mixed_mcmc: design/cluster row mismatch"
      if (size(y_con, 1) /= n .or. size(y_cat, 1) /= n) error stop "jomo1ran_mixed_mcmc: response row mismatch"
      n_con = size(y_con, 2)
      n_cat = size(y_cat, 2)
      if (size(n_levels) /= n_cat) error stop "jomo1ran_mixed_mcmc: n_levels size mismatch"
      if (n_iter <= 0) error stop "jomo1ran_mixed_mcmc: n_iter must be positive"
      if (any(cluster < 1)) error stop "jomo1ran_mixed_mcmc: clusters must be positive one-based labels"
      g_count = maxval(cluster)
      do g = 1, g_count
         if (count(cluster == g) == 0) error stop "jomo1ran_mixed_mcmc: cluster labels must be contiguous"
      end do
      p = latent_dimension(n_con, n_levels)
      qx = size(x, 2)
      qz = size(z, 2)
      d = qz * p
      if (size(prior_scale, 1) /= p .or. size(prior_scale, 2) /= p) error stop "jomo1ran_mixed_mcmc: residual prior shape mismatch"
      if (size(prior_random_scale, 1) /= d .or. size(prior_random_scale, 2) /= d) &
         error stop "jomo1ran_mixed_mcmc: random prior shape mismatch"
      if (n <= max(p, qx)) error stop "jomo1ran_mixed_mcmc: insufficient rows for fixed-effect/covariance update"

      call initialize_latent_data(y_con, con_observed, y_cat, cat_observed, n_levels, result%level1%latent, observed)
      allocate(beta(qx, p), omega(p, p), u(g_count, qz, p), covu(d, d), fixed(p, p))
      allocate(xtx(qx, qx), xtx_inv(qx, qx), bhat(qx, p), random_mean(n, p), mean(n, p), target(n, p))
      allocate(residual(n, p), scale(p, p), covu_precision(d, d), precision(p, p))
      allocate(post_precision(d, d), post_cov(d, d), rhs(d), post_mean(d), uvec(d), hmat(p, d))
      beta = 0.0_dp
      omega = 0.0_dp
      covu = 0.0_dp
      u = 0.0_dp
      do i = 1, p
         omega(i, i) = 1.0_dp
      end do
      do i = 1, d
         covu(i, i) = 1.0_dp
      end do
      if (present(beta_start)) then
         if (any(shape(beta_start) /= shape(beta))) error stop "jomo1ran_mixed_mcmc: beta_start shape mismatch"
         beta = beta_start
      end if
      if (present(omega_start)) then
         if (any(shape(omega_start) /= shape(omega))) error stop "jomo1ran_mixed_mcmc: omega_start shape mismatch"
         omega = omega_start
      end if
      if (present(random_start)) then
         if (any(shape(random_start) /= shape(u))) error stop "jomo1ran_mixed_mcmc: random_start shape mismatch"
         u = random_start
      end if
      if (present(random_cov_start)) then
         if (any(shape(random_cov_start) /= shape(covu))) error stop "jomo1ran_mixed_mcmc: random_cov_start shape mismatch"
         covu = random_cov_start
      end if
      call set_categorical_covariance(omega, n_con, n_levels, fixed)
      if (.not. is_spd(omega)) error stop "jomo1ran_mixed_mcmc: starting omega is not positive definite"
      if (.not. is_spd(covu)) error stop "jomo1ran_mixed_mcmc: starting random covariance is not positive definite"

      xtx = matmul(transpose(x), x)
      call inverse_spd(xtx, xtx_inv, info)
      if (info /= 0) error stop "jomo1ran_mixed_mcmc: fixed-effect design is rank deficient"
      allocate(result%level1%beta_mean(qx, p), result%level1%omega_mean(p, p))
      allocate(result%random_effects_mean(g_count, qz, p), result%random_covariance_mean(d, d))
      result%level1%beta_mean = 0.0_dp
      result%level1%omega_mean = 0.0_dp
      result%random_effects_mean = 0.0_dp
      result%random_covariance_mean = 0.0_dp
      total_failures = 0

      do iter = 1, n_iter
         call compute_random_mean(z, cluster, u, random_mean)
         mean = matmul(x, beta) + random_mean
         if (n_cat > 0) then
            call update_observed_categories(rng, y_cat, cat_observed, n_levels, n_con, mean, omega, &
               result%level1%latent, failures)
            total_failures = total_failures + failures
         end if

         target = result%level1%latent - random_mean
         bhat = matmul(xtx_inv, matmul(transpose(x), target))
         call matrix_normal_sample(rng, bhat, xtx_inv, omega, beta, info)
         if (info /= 0) error stop "jomo1ran_mixed_mcmc: beta draw failed"

         call inverse_spd(omega, precision, info)
         if (info /= 0) error stop "jomo1ran_mixed_mcmc: omega inversion failed"
         call inverse_spd(covu, covu_precision, info)
         if (info /= 0) error stop "jomo1ran_mixed_mcmc: random covariance inversion failed"
         do g = 1, g_count
            post_precision = covu_precision
            rhs = 0.0_dp
            do i = 1, n
               if (cluster(i) /= g) cycle
               call random_design_row(z(i, :), p, hmat)
               post_precision = post_precision + matmul(transpose(hmat), matmul(precision, hmat))
               rhs = rhs + matmul(transpose(hmat), matmul(precision, result%level1%latent(i, :) - matmul(x(i, :), beta)))
            end do
            call inverse_spd(post_precision, post_cov, info)
            if (info /= 0) error stop "jomo1ran_mixed_mcmc: random-effect posterior inversion failed"
            post_mean = matmul(post_cov, rhs)
            call mvnormal_sample(rng, post_mean, post_cov, uvec, info)
            if (info /= 0) error stop "jomo1ran_mixed_mcmc: random-effect draw failed"
            call unpack_random_vector(uvec, u(g, :, :))
         end do

         scale = prior_random_scale
         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            scale = scale + spread(uvec, 2, d) * spread(uvec, 1, d)
         end do
         call invwishart_sample(rng, real(g_count + d, dp), scale, covu, info)
         if (info /= 0) error stop "jomo1ran_mixed_mcmc: random covariance draw failed"

         call compute_random_mean(z, cluster, u, random_mean)
         mean = matmul(x, beta) + random_mean
         residual = result%level1%latent - mean
         if (n_cat == 0) then
            scale = prior_scale + matmul(transpose(residual), residual)
            call invwishart_sample(rng, real(n - 1, dp), scale, omega, info)
            if (info /= 0) error stop "jomo1ran_mixed_mcmc: residual covariance draw failed"
         else
            call mh_covariance_multilevel(rng, residual, fixed, omega)
         end if
         call impute_missing_rows(rng, observed, mean, omega, result%level1%latent, info)
         if (info /= 0) error stop "jomo1ran_mixed_mcmc: conditional imputation failed"

         result%level1%beta_mean = result%level1%beta_mean + beta
         result%level1%omega_mean = result%level1%omega_mean + omega
         result%random_effects_mean = result%random_effects_mean + u
         result%random_covariance_mean = result%random_covariance_mean + covu
      end do

      result%level1%beta_mean = result%level1%beta_mean / real(n_iter, dp)
      result%level1%omega_mean = result%level1%omega_mean / real(n_iter, dp)
      result%random_effects_mean = result%random_effects_mean / real(n_iter, dp)
      result%random_covariance_mean = result%random_covariance_mean / real(n_iter, dp)
      allocate(result%level1%beta(qx, p), result%level1%omega(p, p))
      allocate(result%random_effects(g_count, qz, p), result%random_covariance(d, d))
      result%level1%beta = beta
      result%level1%omega = omega
      result%random_effects = u
      result%random_covariance = covu
      allocate(result%level1%continuous(n, n_con), result%level1%categorical(n, n_cat))
      if (n_con > 0) result%level1%continuous = result%level1%latent(:, 1:n_con)
      if (n_cat > 0) call decode_categories(result%level1%latent, n_con, n_levels, result%level1%categorical)
      result%level1%latent_rejection_failures = total_failures
      result%level1%iterations = n_iter
      do j = 1, n_con
         do i = 1, n
            if (con_observed(i, j)) result%level1%continuous(i, j) = y_con(i, j)
         end do
      end do
      do j = 1, n_cat
         do i = 1, n
            if (cat_observed(i, j)) result%level1%categorical(i, j) = y_cat(i, j)
         end do
      end do
   end subroutine jomo1ran_mixed_mcmc

   subroutine jomo1rancon_mcmc(rng, y, observed, x, z, cluster, n_iter, prior_scale, prior_random_scale, result)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all multilevel MCMC draws.
      real(dp), intent(in) :: y(:, :) !! Continuous level-1 response matrix, shape n by p.
      logical, intent(in) :: observed(:, :) !! Observation mask for y.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by qx.
      real(dp), intent(in) :: z(:, :) !! Fully observed random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels for rows of y.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale(:, :) !! Residual inverse-Wishart scale matrix, p by p.
      real(dp), intent(in) :: prior_random_scale(:, :) !! Random-effect inverse-Wishart scale matrix, qz*p by qz*p.
      type(jomo1ran_result), intent(out) :: result !! Final imputation and multilevel parameter draws/posterior means.
      integer, allocatable :: y_cat(:, :)
      integer, allocatable :: n_levels(:)
      logical, allocatable :: cat_observed(:, :)

      allocate(y_cat(size(y, 1), 0), cat_observed(size(y, 1), 0), n_levels(0))
      call jomo1ran_mixed_mcmc(rng, y, observed, y_cat, cat_observed, n_levels, x, z, cluster, n_iter, prior_scale, &
         prior_random_scale, result)
   end subroutine jomo1rancon_mcmc

   pure subroutine compute_random_mean(z, cluster, u, random_mean)
      real(dp), intent(in) :: z(:, :) !! Random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based cluster label for each design row.
      real(dp), intent(in) :: u(:, :, :) !! Random effects, shape G by qz by p.
      real(dp), intent(out) :: random_mean(:, :) !! Row-specific random-effect contribution, shape n by p.
      integer :: i
      integer :: j

      if (size(random_mean, 1) /= size(z, 1) .or. size(random_mean, 2) /= size(u, 3)) &
         error stop "compute_random_mean: shape mismatch"
      do i = 1, size(z, 1)
         do j = 1, size(u, 3)
            random_mean(i, j) = dot_product(z(i, :), u(cluster(i), :, j))
         end do
      end do
   end subroutine compute_random_mean

   pure subroutine random_design_row(zrow, p, hmat)
      real(dp), intent(in) :: zrow(:) !! One random-effect design row of length qz.
      integer, intent(in) :: p !! Number of latent response components.
      real(dp), intent(out) :: hmat(:, :) !! Block-diagonal observation design, shape p by qz*p, for vec(random effects).
      integer :: j
      integer :: qz

      qz = size(zrow)
      if (size(hmat, 1) /= p .or. size(hmat, 2) /= qz * p) error stop "random_design_row: shape mismatch"
      hmat = 0.0_dp
      do j = 1, p
         hmat(j, (j - 1) * qz + 1:j * qz) = zrow
      end do
   end subroutine random_design_row

   pure subroutine pack_random_matrix(u, v)
      real(dp), intent(in) :: u(:, :) !! Random-effect matrix, shape qz by p.
      real(dp), intent(out) :: v(:) !! Response-major vectorization of u, length qz*p.
      integer :: j
      integer :: qz

      qz = size(u, 1)
      if (size(v) /= size(u)) error stop "pack_random_matrix: shape mismatch"
      do j = 1, size(u, 2)
         v((j - 1) * qz + 1:j * qz) = u(:, j)
      end do
   end subroutine pack_random_matrix

   pure subroutine unpack_random_vector(v, u)
      real(dp), intent(in) :: v(:) !! Response-major vectorized random effects, length qz*p.
      real(dp), intent(out) :: u(:, :) !! Random-effect matrix, shape qz by p.
      integer :: j
      integer :: qz

      qz = size(u, 1)
      if (size(v) /= size(u)) error stop "unpack_random_vector: shape mismatch"
      do j = 1, size(u, 2)
         u(:, j) = v((j - 1) * qz + 1:j * qz)
      end do
   end subroutine unpack_random_vector

   subroutine mh_covariance_multilevel(rng, residual, fixed, omega)
      use jomo_rng, only : rng_uniform, rng_normal
      use jomo_linalg, only : logdet_spd, quadratic_form
      type(rng_state), intent(inout) :: rng !! Mutable generator state for covariance proposals and acceptance draws.
      real(dp), intent(in) :: residual(:, :) !! Current complete-data residual matrix, shape n by p.
      logical, intent(in) :: fixed(:, :) !! Identification mask for fixed categorical covariance elements.
      real(dp), intent(inout) :: omega(:, :) !! Residual covariance matrix updated in-place.
      integer :: n
      integer :: p
      integer :: j
      integer :: k
      integer :: tries
      integer :: info
      integer :: i
      real(dp) :: current_ll
      real(dp) :: proposed_ll
      real(dp) :: sd
      real(dp) :: candidate
      real(dp), allocatable :: proposal(:, :)
      real(dp), allocatable :: inv(:, :)
      logical :: found

      n = size(residual, 1)
      p = size(residual, 2)
      allocate(proposal(p, p), inv(p, p))
      call inverse_spd(omega, inv, info)
      if (info /= 0) error stop "mh_covariance_multilevel: omega inversion failed"
      current_ll = -0.5_dp * real(n, dp) * logdet_spd(omega)
      do i = 1, n
         current_ll = current_ll - 0.5_dp * quadratic_form(residual(i, :), inv)
      end do
      do k = 1, p
         do j = 1, k
            if (fixed(j, k)) cycle
            if (j == k) then
               sd = omega(j, j) * sqrt(11.6_dp / real(n, dp))
            else
               sd = 0.1_dp * sqrt(omega(j, j) * omega(k, k))
            end if
            if (sd <= 0.0_dp) cycle
            found = .false.
            do tries = 1, 100
               candidate = rng_normal(rng, omega(j, k), sd)
               proposal = omega
               proposal(j, k) = candidate
               proposal(k, j) = candidate
               if (is_spd(proposal)) then
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) cycle
            call inverse_spd(proposal, inv, info)
            if (info /= 0) cycle
            proposed_ll = -0.5_dp * real(n, dp) * logdet_spd(proposal)
            do i = 1, n
               proposed_ll = proposed_ll - 0.5_dp * quadratic_form(residual(i, :), inv)
            end do
            if (log(rng_uniform(rng)) < proposed_ll - current_ll) then
               omega = proposal
               current_ll = proposed_ll
            end if
         end do
      end do
   end subroutine mh_covariance_multilevel

end module jomo_multilevel
