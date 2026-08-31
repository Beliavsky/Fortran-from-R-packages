! Heteroscedastic multilevel joint-model MCMC kernels for jomo.
! This module translates the main numerical structure of jomo1ranconhr,
! jomo1rancathr and jomo1ranmixhr using cluster-specific residual covariance.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_heteroscedastic
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state, rng_uniform, rng_normal, rng_student_t
   use jomo_linalg, only : inverse_spd, is_spd, logdet_spd, quadratic_form
   use jomo_distributions, only : mvnormal_sample, invwishart_sample, wishart_sample
   use jomo_latent, only : latent_dimension, initialize_latent_data, set_categorical_covariance
   use jomo_latent, only : decode_categories, update_observed_categories, impute_missing_rows
   use jomo_single_level, only : jomo1_result
   use jomo_multilevel, only : compute_random_mean, random_design_row, pack_random_matrix, unpack_random_vector
   implicit none
   private

   type, public :: jomo1ranhr_result
      type(jomo1_result) :: level1
      real(dp), allocatable :: random_effects(:, :, :)
      real(dp), allocatable :: random_covariance(:, :)
      real(dp), allocatable :: cluster_covariance(:, :, :)
      real(dp), allocatable :: hierarchy_scale(:, :)
      real(dp) :: hierarchy_df = 0.0_dp
   end type jomo1ranhr_result

   public :: jomo1ranhr_mixed_mcmc
   public :: jomo1ranconhr_mcmc
   public :: update_hierarchy_degrees
   public :: mh_cluster_covariance

contains

   subroutine jomo1ranhr_mixed_mcmc(rng, y_con, con_observed, y_cat, cat_observed, n_levels, x, z, cluster, n_iter, &
      prior_scale, prior_random_scale, hierarchy_df, result, hierarchical, update_hierarchy_df, hierarchy_df_prior)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all heterogeneous-model MCMC draws.
      real(dp), intent(in) :: y_con(:, :) !! Continuous level-1 responses, shape n by n_con; masked entries are ignored.
      logical, intent(in) :: con_observed(:, :) !! Observation mask for continuous level-1 responses.
      integer, intent(in) :: y_cat(:, :) !! Categorical level-1 responses encoded 1..K, shape n by n_cat.
      logical, intent(in) :: cat_observed(:, :) !! Observation mask for categorical level-1 responses.
      integer, intent(in) :: n_levels(:) !! Number of levels K for each categorical level-1 response.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by qx.
      real(dp), intent(in) :: z(:, :) !! Fully observed random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels for level-1 rows.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale(:, :) !! Positive-definite p by p residual/hierarchy scale matrix (upstream Sp).
      real(dp), intent(in) :: prior_random_scale(:, :) !! Positive-definite qz*p by qz*p inverse-Wishart scale for random effects.
      real(dp), intent(in) :: hierarchy_df !! Initial inverse-Wishart degrees parameter a for cluster-specific covariance matrices.
      type(jomo1ranhr_result), intent(out) :: result !! Final heterogeneous imputation and parameter state.
      logical, intent(in), optional :: hierarchical !! True for random-covariance hierarchy; false for fixed cluster covariances.
      logical, intent(in), optional :: update_hierarchy_df !! If true, sample a by the upstream Newton-centered Student-t MH step.
      real(dp), intent(in), optional :: hierarchy_df_prior !! Chi-square prior degrees eta for a; default is latent dimension p.
      integer :: n
      integer :: n_con
      integer :: n_cat
      integer :: p
      integer :: qx
      integer :: qz
      integer :: d
      integer :: g_count
      integer :: iter
      integer :: i
      integer :: j
      integer :: g
      integer :: info
      integer :: failures
      integer :: n_g
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: beta_vec(:)
      real(dp), allocatable :: omega(:, :, :)
      real(dp), allocatable :: u(:, :, :)
      real(dp), allocatable :: covu(:, :)
      real(dp), allocatable :: covu_precision(:, :)
      real(dp), allocatable :: random_mean(:, :)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: hbeta(:, :)
      real(dp), allocatable :: hu(:, :)
      real(dp), allocatable :: post_precision(:, :)
      real(dp), allocatable :: post_cov(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: post_mean(:)
      real(dp), allocatable :: uvec(:)
      real(dp), allocatable :: scale(:, :)
      real(dp), allocatable :: sum_precision(:, :)
      real(dp), allocatable :: prior_precision(:, :)
      real(dp), allocatable :: hierarchy_post_scale(:, :)
      real(dp), allocatable :: hierarchy_scale(:, :)
      logical, allocatable :: observed(:, :)
      logical, allocatable :: fixed(:, :)
      logical :: use_hierarchy
      logical :: sample_hierarchy_df
      real(dp) :: a_current
      real(dp) :: df_prior

      n = size(x, 1)
      if (size(z, 1) /= n .or. size(cluster) /= n) error stop "jomo1ranhr_mixed_mcmc: design/cluster row mismatch"
      if (size(y_con, 1) /= n .or. size(y_cat, 1) /= n) error stop "jomo1ranhr_mixed_mcmc: response row mismatch"
      n_con = size(y_con, 2)
      n_cat = size(y_cat, 2)
      if (size(n_levels) /= n_cat) error stop "jomo1ranhr_mixed_mcmc: n_levels size mismatch"
      if (n_iter <= 0) error stop "jomo1ranhr_mixed_mcmc: n_iter must be positive"
      if (hierarchy_df <= 0.0_dp) error stop "jomo1ranhr_mixed_mcmc: hierarchy_df must be positive"
      if (any(cluster < 1)) error stop "jomo1ranhr_mixed_mcmc: cluster labels must be positive"
      g_count = maxval(cluster)
      do g = 1, g_count
         if (count(cluster == g) == 0) error stop "jomo1ranhr_mixed_mcmc: cluster labels must be contiguous"
      end do
      p = latent_dimension(n_con, n_levels)
      qx = size(x, 2)
      qz = size(z, 2)
      d = qz * p
      if (size(prior_scale, 1) /= p .or. size(prior_scale, 2) /= p) error stop "jomo1ranhr_mixed_mcmc: prior scale shape mismatch"
      if (size(prior_random_scale, 1) /= d .or. size(prior_random_scale, 2) /= d) &
         error stop "jomo1ranhr_mixed_mcmc: random prior shape mismatch"
      use_hierarchy = .true.
      if (present(hierarchical)) use_hierarchy = hierarchical
      sample_hierarchy_df = use_hierarchy
      if (present(update_hierarchy_df)) sample_hierarchy_df = update_hierarchy_df .and. use_hierarchy
      a_current = hierarchy_df
      df_prior = real(p, dp)
      if (present(hierarchy_df_prior)) df_prior = hierarchy_df_prior
      if (df_prior <= 0.0_dp) error stop "jomo1ranhr_mixed_mcmc: hierarchy_df_prior must be positive"
      if (sample_hierarchy_df .and. a_current < real(p, dp)) &
         error stop "jomo1ranhr_mixed_mcmc: sampled hierarchy_df must be at least p"

      call initialize_latent_data(y_con, con_observed, y_cat, cat_observed, n_levels, result%level1%latent, observed)
      allocate(beta(qx, p), beta_vec(qx * p), omega(g_count, p, p), u(g_count, qz, p), covu(d, d))
      allocate(covu_precision(d, d), random_mean(n, p), mean(n, p), residual(n, p), precision(p, p))
      allocate(hbeta(p, qx * p), hu(p, d), post_precision(max(qx * p, d), max(qx * p, d)))
      allocate(post_cov(max(qx * p, d), max(qx * p, d)), rhs(max(qx * p, d)), post_mean(max(qx * p, d)))
      allocate(uvec(d), fixed(p, p), sum_precision(p, p), prior_precision(p, p), hierarchy_post_scale(p, p))
      allocate(hierarchy_scale(p, p))
      beta = 0.0_dp
      u = 0.0_dp
      covu = 0.0_dp
      do i = 1, d
         covu(i, i) = 1.0_dp
      end do
      do g = 1, g_count
         omega(g, :, :) = 0.0_dp
         do i = 1, p
            omega(g, i, i) = 1.0_dp
         end do
         call set_categorical_covariance(omega(g, :, :), n_con, n_levels, fixed)
      end do
      hierarchy_scale = prior_scale
      call inverse_spd(prior_scale, prior_precision, info)
      if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: prior_scale is not positive definite"
      if (.not. is_spd(covu)) error stop "jomo1ranhr_mixed_mcmc: initial random covariance failure"

      do iter = 1, n_iter
         call compute_random_mean(z, cluster, u, random_mean)
         mean = matmul(x, beta) + random_mean
         if (n_cat > 0) then
            do i = 1, n
               call update_observed_categories(rng, y_cat(i:i, :), cat_observed(i:i, :), n_levels, n_con, &
                  mean(i:i, :), omega(cluster(i), :, :), result%level1%latent(i:i, :), failures)
            end do
         end if

         post_precision(1:qx * p, 1:qx * p) = 0.0_dp
         rhs(1:qx * p) = 0.0_dp
         do i = 1, n
            call inverse_spd(omega(cluster(i), :, :), precision, info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: cluster covariance inversion failed"
            call fixed_design_row(x(i, :), p, hbeta)
            post_precision(1:qx * p, 1:qx * p) = post_precision(1:qx * p, 1:qx * p) + &
               matmul(transpose(hbeta), matmul(precision, hbeta))
            rhs(1:qx * p) = rhs(1:qx * p) + matmul(transpose(hbeta), &
               matmul(precision, result%level1%latent(i, :) - random_mean(i, :)))
         end do
         call inverse_spd(post_precision(1:qx * p, 1:qx * p), post_cov(1:qx * p, 1:qx * p), info)
         if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: beta posterior inversion failed"
         post_mean(1:qx * p) = matmul(post_cov(1:qx * p, 1:qx * p), rhs(1:qx * p))
         call mvnormal_sample(rng, post_mean(1:qx * p), post_cov(1:qx * p, 1:qx * p), beta_vec, info)
         if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: beta draw failed"
         call unpack_response_major(beta_vec, beta)

         call inverse_spd(covu, covu_precision, info)
         if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: random covariance inversion failed"
         do g = 1, g_count
            call inverse_spd(omega(g, :, :), precision, info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: cluster precision failure"
            post_precision(1:d, 1:d) = covu_precision
            rhs(1:d) = 0.0_dp
            do i = 1, n
               if (cluster(i) /= g) cycle
               call random_design_row(z(i, :), p, hu)
               post_precision(1:d, 1:d) = post_precision(1:d, 1:d) + matmul(transpose(hu), matmul(precision, hu))
               rhs(1:d) = rhs(1:d) + matmul(transpose(hu), &
                  matmul(precision, result%level1%latent(i, :) - matmul(x(i, :), beta)))
            end do
            call inverse_spd(post_precision(1:d, 1:d), post_cov(1:d, 1:d), info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: random-effect posterior inversion failed"
            post_mean(1:d) = matmul(post_cov(1:d, 1:d), rhs(1:d))
            call mvnormal_sample(rng, post_mean(1:d), post_cov(1:d, 1:d), uvec, info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: random-effect draw failed"
            call unpack_random_vector(uvec, u(g, :, :))
         end do

         scale = prior_random_scale
         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            scale = scale + spread(uvec, 2, d) * spread(uvec, 1, d)
         end do
         call invwishart_sample(rng, real(g_count + d, dp), scale, covu, info)
         if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: random covariance draw failed"

         call compute_random_mean(z, cluster, u, random_mean)
         mean = matmul(x, beta) + random_mean
         residual = result%level1%latent - mean

         if (use_hierarchy) then
            sum_precision = prior_precision
            do g = 1, g_count
               call inverse_spd(omega(g, :, :), precision, info)
               if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: hierarchy precision inversion failed"
               sum_precision = sum_precision + precision
            end do
            call inverse_spd(sum_precision, hierarchy_post_scale, info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: hierarchy scale inversion failed"
            call wishart_sample(rng, real(g_count, dp) * a_current + real(p + 1, dp), hierarchy_post_scale, &
               hierarchy_scale, info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: hierarchy scale draw failed"
            if (sample_hierarchy_df) call update_hierarchy_degrees(rng, a_current, df_prior, omega, sum_precision)
         else
            hierarchy_scale = prior_scale
         end if

         do g = 1, g_count
            n_g = count(cluster == g)
            scale = 0.0_dp * hierarchy_scale
            do i = 1, n
               if (cluster(i) /= g) cycle
               scale = scale + spread(residual(i, :), 2, p) * spread(residual(i, :), 1, p)
            end do
            if (n_cat == 0) then
               if (use_hierarchy) then
                  scale = scale + hierarchy_scale
                  call invwishart_sample(rng, real(n_g, dp) + a_current, scale, omega(g, :, :), info)
               else
                  scale = scale + prior_scale
                  call invwishart_sample(rng, real(n_g, dp), scale, omega(g, :, :), info)
               end if
               if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: cluster covariance draw failed"
            else
               call mh_cluster_covariance(rng, residual, cluster, g, fixed, hierarchy_scale, a_current, use_hierarchy, &
                  omega(g, :, :))
            end if
         end do

         do i = 1, n
            call impute_missing_rows(rng, observed(i:i, :), mean(i:i, :), omega(cluster(i), :, :), &
               result%level1%latent(i:i, :), info)
            if (info /= 0) error stop "jomo1ranhr_mixed_mcmc: conditional imputation failed"
         end do
      end do

      allocate(result%level1%beta(qx, p), result%level1%omega(p, p))
      allocate(result%random_effects(g_count, qz, p), result%random_covariance(d, d))
      allocate(result%cluster_covariance(g_count, p, p), result%hierarchy_scale(p, p))
      result%level1%beta = beta
      result%level1%omega = sum(omega, dim=1) / real(g_count, dp)
      result%random_effects = u
      result%random_covariance = covu
      result%cluster_covariance = omega
      result%hierarchy_scale = hierarchy_scale
      result%hierarchy_df = a_current
      allocate(result%level1%continuous(n, n_con), result%level1%categorical(n, n_cat))
      if (n_con > 0) result%level1%continuous = result%level1%latent(:, 1:n_con)
      if (n_cat > 0) call decode_categories(result%level1%latent, n_con, n_levels, result%level1%categorical)
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
   end subroutine jomo1ranhr_mixed_mcmc

   subroutine jomo1ranconhr_mcmc(rng, y, observed, x, z, cluster, n_iter, prior_scale, prior_random_scale, hierarchy_df, &
      result, update_hierarchy_df, hierarchy_df_prior)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all heterogeneous continuous-model draws.
      real(dp), intent(in) :: y(:, :) !! Continuous level-1 response matrix, shape n by p.
      logical, intent(in) :: observed(:, :) !! Observation mask for y.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by qx.
      real(dp), intent(in) :: z(:, :) !! Fully observed random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale(:, :) !! Positive-definite p by p residual/hierarchy scale matrix.
      real(dp), intent(in) :: prior_random_scale(:, :) !! Positive-definite qz*p by qz*p random-effect scale matrix.
      real(dp), intent(in) :: hierarchy_df !! Initial inverse-Wishart degrees parameter a for cluster-specific covariance matrices.
      type(jomo1ranhr_result), intent(out) :: result !! Final heterogeneous continuous imputation and parameter state.
      logical, intent(in), optional :: update_hierarchy_df !! If true, sample a with the upstream random-covariance MH update.
      real(dp), intent(in), optional :: hierarchy_df_prior !! Chi-square prior degrees eta for a; default is response dimension.
      integer, allocatable :: y_cat(:, :)
      integer, allocatable :: n_levels(:)
      logical, allocatable :: cat_observed(:, :)

      allocate(y_cat(size(y, 1), 0), cat_observed(size(y, 1), 0), n_levels(0))
      call jomo1ranhr_mixed_mcmc(rng, y, observed, y_cat, cat_observed, n_levels, x, z, cluster, n_iter, prior_scale, &
         prior_random_scale, hierarchy_df, result, update_hierarchy_df=update_hierarchy_df, hierarchy_df_prior=hierarchy_df_prior)
   end subroutine jomo1ranconhr_mcmc

   pure subroutine fixed_design_row(xrow, p, hmat)
      real(dp), intent(in) :: xrow(:) !! One fixed-effect design row of length qx.
      integer, intent(in) :: p !! Number of latent response components.
      real(dp), intent(out) :: hmat(:, :) !! Block-diagonal design mapping response-major vec(beta) to a p-vector.
      integer :: j
      integer :: qx

      qx = size(xrow)
      if (size(hmat, 1) /= p .or. size(hmat, 2) /= qx * p) error stop "fixed_design_row: shape mismatch"
      hmat = 0.0_dp
      do j = 1, p
         hmat(j, (j - 1) * qx + 1:j * qx) = xrow
      end do
   end subroutine fixed_design_row

   pure subroutine unpack_response_major(v, b)
      real(dp), intent(in) :: v(:) !! Response-major vectorization of a q by p coefficient matrix.
      real(dp), intent(out) :: b(:, :) !! Coefficient matrix with q rows and p response columns.
      integer :: j
      integer :: q

      q = size(b, 1)
      if (size(v) /= size(b)) error stop "unpack_response_major: shape mismatch"
      do j = 1, size(b, 2)
         b(:, j) = v((j - 1) * q + 1:j * q)
      end do
   end subroutine unpack_response_major

   subroutine update_hierarchy_degrees(rng, a, eta, omega, sum_precision)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the scalar degrees-of-freedom Metropolis proposal.
      real(dp), intent(inout) :: a !! Current inverse-Wishart degrees parameter; updated in place and constrained to a >= p.
      real(dp), intent(in) :: eta !! Chi-square prior degrees for a, matching upstream a.prior.
      real(dp), intent(in) :: omega(:, :, :) !! Current cluster covariance matrices, shape n_cluster by p by p.
      real(dp), intent(in) :: sum_precision(:, :) !! Prior precision plus the sum of current cluster precision matrices.
      integer, parameter :: max_newton = 50
      real(dp), parameter :: dx = 0.001_dp
      real(dp), parameter :: precision = 0.001_dp
      integer :: iter
      integer :: p
      real(dp) :: u_current
      real(dp) :: u_mode
      real(dp) :: u_next
      real(dp) :: u_prop
      real(dp) :: deriv1
      real(dp) :: deriv2
      real(dp) :: lambda
      real(dp) :: log_alpha
      real(dp) :: a_prop
      logical :: converged

      p = size(omega, 2)
      u_current = log(a + real(p, dp))
      u_mode = u_current
      converged = .false.
      do iter = 1, max_newton
         call hierarchy_target_derivatives(eta, u_mode, omega, sum_precision, dx, deriv1, deriv2)
         if (.not. ieee_is_finite(deriv1) .or. .not. ieee_is_finite(deriv2)) exit
         if (abs(deriv2) <= tiny(1.0_dp)) exit
         u_next = u_mode - deriv1 / deriv2
         if (.not. ieee_is_finite(u_next)) exit
         if (abs(u_next - u_mode) <= precision * max(1.0_dp, abs(u_next))) then
            u_mode = u_next
            converged = .true.
            exit
         end if
         u_mode = u_next
      end do
      if (.not. converged) u_mode = u_current

      call hierarchy_target_derivatives(eta, u_mode, omega, sum_precision, dx, deriv1, deriv2)
      if (.not. ieee_is_finite(deriv2) .or. deriv2 >= 0.0_dp) return
      lambda = sqrt(-5.0_dp / (4.0_dp * deriv2))
      if (.not. ieee_is_finite(lambda) .or. lambda <= 0.0_dp) return
      u_prop = u_mode + lambda * rng_student_t(rng, 4.0_dp)
      if (.not. ieee_is_finite(u_prop)) return
      a_prop = exp(u_prop) - real(p, dp)
      if (.not. ieee_is_finite(a_prop) .or. a_prop < real(p, dp)) return

      log_alpha = hierarchy_df_logtarget(eta, u_prop, omega, sum_precision) - &
         hierarchy_df_logtarget(eta, u_current, omega, sum_precision)
      log_alpha = log_alpha + log_student_proposal_kernel(u_current, u_mode, lambda) - &
         log_student_proposal_kernel(u_prop, u_mode, lambda)
      if (log(rng_uniform(rng)) < min(0.0_dp, log_alpha)) a = a_prop
   end subroutine update_hierarchy_degrees

   subroutine hierarchy_target_derivatives(eta, u, omega, sum_precision, dx, first, second)
      real(dp), intent(in) :: eta !! Chi-square prior degrees for the inverse-Wishart degrees parameter.
      real(dp), intent(in) :: u !! Transformed parameter u = log(a + p) where derivatives are evaluated.
      real(dp), intent(in) :: omega(:, :, :) !! Current cluster covariance matrices.
      real(dp), intent(in) :: sum_precision(:, :) !! Prior precision plus current cluster precision sum.
      real(dp), intent(in) :: dx !! Positive central finite-difference increment in transformed-parameter units.
      real(dp), intent(out) :: first !! Central first derivative of the transformed log target.
      real(dp), intent(out) :: second !! Central second derivative of the transformed log target.
      real(dp) :: fm
      real(dp) :: f0
      real(dp) :: fp

      fm = hierarchy_df_logtarget(eta, u - dx, omega, sum_precision)
      f0 = hierarchy_df_logtarget(eta, u, omega, sum_precision)
      fp = hierarchy_df_logtarget(eta, u + dx, omega, sum_precision)
      first = (fp - fm) / (2.0_dp * dx)
      second = (fp - 2.0_dp * f0 + fm) / (dx * dx)
   end subroutine hierarchy_target_derivatives

   function hierarchy_df_logtarget(eta, u, omega, sum_precision) result(value)
      real(dp), intent(in) :: eta !! Chi-square prior degrees for a.
      real(dp), intent(in) :: u !! Transformed degrees parameter u = log(a + p).
      real(dp), intent(in) :: omega(:, :, :) !! Current cluster covariance matrices.
      real(dp), intent(in) :: sum_precision(:, :) !! Prior precision plus current cluster precision sum.
      real(dp) :: value
      integer :: g
      integer :: g_count
      integer :: p
      real(dp) :: a
      real(dp) :: logdet_precision
      real(dp) :: logdet_sum_precision

      p = size(omega, 2)
      g_count = size(omega, 1)
      a = exp(u) - real(p, dp)
      if (.not. ieee_is_finite(a) .or. a <= real(p - 1, dp)) then
         value = -huge(1.0_dp)
         return
      end if
      logdet_sum_precision = logdet_spd(sum_precision)
      if (.not. ieee_is_finite(logdet_sum_precision)) then
         value = -huge(1.0_dp)
         return
      end if

      value = log_chisq_density(eta, a)
      value = value - real(g_count, dp) * log_multivariate_gamma_kernel(p, 0.5_dp * a)
      do g = 1, g_count
         logdet_precision = -logdet_spd(omega(g, :, :))
         if (.not. ieee_is_finite(logdet_precision)) then
            value = -huge(1.0_dp)
            return
         end if
         value = value + 0.5_dp * (a + real(p + 1, dp)) * logdet_precision
      end do
      value = value - 0.5_dp * (a * real(g_count, dp) + real(p + 1, dp)) * logdet_sum_precision
      value = value + log_multivariate_gamma_kernel(p, &
         0.5_dp * (a * real(g_count, dp) + real(p + 1, dp))) + u
   end function hierarchy_df_logtarget

   pure function log_chisq_density(df, x) result(value)
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
      real(dp), intent(in) :: x !! Positive chi-square variate at which the log density is evaluated.
      real(dp) :: value

      if (df <= 0.0_dp .or. x <= 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = (0.5_dp * df - 1.0_dp) * log(x) - 0.5_dp * x - &
            0.5_dp * df * log(2.0_dp) - log_gamma(0.5_dp * df)
      end if
   end function log_chisq_density

   pure function log_multivariate_gamma_kernel(p, a) result(value)
      integer, intent(in) :: p !! Positive multivariate-gamma dimension.
      real(dp), intent(in) :: a !! Multivariate-gamma shape; must exceed (p-1)/2.
      real(dp) :: value
      integer :: j

      if (p <= 0 .or. a <= 0.5_dp * real(p - 1, dp)) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do j = 1, p
         value = value + log_gamma(a + 0.5_dp * real(1 - j, dp))
      end do
   end function log_multivariate_gamma_kernel

   pure function log_student_proposal_kernel(u, u_mode, lambda) result(value)
      real(dp), intent(in) :: u !! Transformed degrees parameter at which the proposal kernel is evaluated.
      real(dp), intent(in) :: u_mode !! Newton mode used as the center of the four-degree Student-t proposal.
      real(dp), intent(in) :: lambda !! Positive proposal scale.
      real(dp) :: value

      value = -2.5_dp * log(1.0_dp + ((u - u_mode) * (u - u_mode)) / (4.0_dp * lambda * lambda))
   end function log_student_proposal_kernel

   subroutine mh_cluster_covariance(rng, residual, cluster, group, fixed, hierarchy_scale, hierarchy_df, use_hierarchy, omega)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for cluster covariance Metropolis updates.
      real(dp), intent(in) :: residual(:, :) !! Complete-data residual matrix for all rows.
      integer, intent(in) :: cluster(:) !! One-based cluster labels associated with residual rows.
      integer, intent(in) :: group !! Cluster whose covariance matrix is being updated.
      logical, intent(in) :: fixed(:, :) !! Identification mask for latent-category covariance elements.
      real(dp), intent(in) :: hierarchy_scale(:, :) !! Current inverse-Wishart scale matrix used when use_hierarchy is true.
      real(dp), intent(in) :: hierarchy_df !! Current inverse-Wishart degrees parameter a used when use_hierarchy is true.
      logical, intent(in) :: use_hierarchy !! Whether the hierarchical inverse-Wishart term enters the target density.
      real(dp), intent(inout) :: omega(:, :) !! Cluster covariance matrix updated in-place.
      integer :: n_g
      integer :: p
      integer :: j
      integer :: k
      integer :: tries
      real(dp) :: current_ll
      real(dp) :: proposed_ll
      real(dp) :: sd
      real(dp) :: candidate
      real(dp), allocatable :: proposal(:, :)
      logical :: found

      n_g = count(cluster == group)
      p = size(omega, 1)
      allocate(proposal(p, p))
      current_ll = cluster_logtarget(residual, cluster, group, omega, hierarchy_scale, hierarchy_df, use_hierarchy)
      do k = 1, p
         do j = 1, k
            if (fixed(j, k)) cycle
            if (j == k) then
               sd = omega(j, j) * sqrt(11.6_dp / real(max(1, n_g), dp))
            else
               sd = 0.1_dp * sqrt(omega(j, j) * omega(k, k))
            end if
            if (sd <= 0.0_dp) cycle
            found = .false.
            do tries = 1, 100
               proposal = omega
               candidate = rng_normal(rng, omega(j, k), sd)
               proposal(j, k) = candidate
               proposal(k, j) = candidate
               if (is_spd(proposal)) then
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) cycle
            proposed_ll = cluster_logtarget(residual, cluster, group, proposal, hierarchy_scale, hierarchy_df, use_hierarchy)
            if (log(rng_uniform(rng)) < proposed_ll - current_ll) then
               omega = proposal
               current_ll = proposed_ll
            end if
         end do
      end do
   end subroutine mh_cluster_covariance

   function cluster_logtarget(residual, cluster, group, omega, hierarchy_scale, hierarchy_df, use_hierarchy) result(value)
      real(dp), intent(in) :: residual(:, :) !! Complete-data residual matrix for all observations.
      integer, intent(in) :: cluster(:) !! One-based cluster labels associated with residual rows.
      integer, intent(in) :: group !! Cluster whose covariance log target is evaluated.
      real(dp), intent(in) :: omega(:, :) !! Candidate cluster covariance matrix.
      real(dp), intent(in) :: hierarchy_scale(:, :) !! Hierarchical inverse-Wishart scale matrix.
      real(dp), intent(in) :: hierarchy_df !! Hierarchical inverse-Wishart degrees parameter a.
      logical, intent(in) :: use_hierarchy !! Whether the inverse-Wishart hierarchy contributes to the target density.
      real(dp) :: value
      real(dp), allocatable :: inv(:, :)
      real(dp), allocatable :: scatter(:, :)
      integer :: i
      integer :: info
      integer :: n_g
      integer :: p

      p = size(omega, 1)
      n_g = count(cluster == group)
      allocate(inv(p, p), scatter(p, p))
      call inverse_spd(omega, inv, info)
      if (info /= 0) then
         value = -huge(1.0_dp)
         return
      end if
      scatter = 0.0_dp
      do i = 1, size(residual, 1)
         if (cluster(i) /= group) cycle
         scatter = scatter + spread(residual(i, :), 2, p) * spread(residual(i, :), 1, p)
      end do
      if (use_hierarchy) then
         scatter = scatter + hierarchy_scale
         value = -0.5_dp * sum(scatter * transpose(inv))
         value = value - 0.5_dp * (real(n_g + p + 1, dp) + hierarchy_df) * logdet_spd(omega)
      else
         value = -0.5_dp * sum(scatter * transpose(inv)) - 0.5_dp * real(n_g, dp) * logdet_spd(omega)
      end if
   end function cluster_logtarget

end module jomo_heteroscedastic
