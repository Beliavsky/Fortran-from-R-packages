! Single-level joint-model multiple-imputation MCMC kernels for jomo.
! This module covers the computational core shared by jomo1con, jomo1cat,
! jomo1mix and their MCMC-chain variants without R data-frame/S3 plumbing.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_single_level
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state, rng_uniform, rng_normal
   use jomo_linalg, only : inverse_spd, logdet_spd, quadratic_form, is_spd
   use jomo_distributions, only : matrix_normal_sample, invwishart_sample
   use jomo_latent, only : latent_dimension, initialize_latent_data, set_categorical_covariance
   use jomo_latent, only : decode_categories, update_observed_categories, impute_missing_rows
   implicit none
   private

   type, public :: jomo1_result
      real(dp), allocatable :: continuous(:, :)
      integer, allocatable :: categorical(:, :)
      real(dp), allocatable :: latent(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: omega(:, :)
      real(dp), allocatable :: beta_mean(:, :)
      real(dp), allocatable :: omega_mean(:, :)
      real(dp), allocatable :: beta_chain(:, :, :)
      real(dp), allocatable :: omega_chain(:, :, :)
      integer :: latent_rejection_failures = 0
      integer :: iterations = 0
   end type jomo1_result

   public :: jomo1_mixed_mcmc
   public :: jomo1con_mcmc
   public :: jomo1cat_mcmc

contains

   subroutine jomo1_mixed_mcmc(rng, y_con, con_observed, y_cat, cat_observed, n_levels, x, n_iter, prior_scale, &
      result, beta_start, omega_start, store_chain)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all Gibbs and Metropolis draws.
      real(dp), intent(in) :: y_con(:, :) !! Continuous responses, shape n by n_con; masked entries may contain any value.
      logical, intent(in) :: con_observed(:, :) !! Observation mask for continuous responses, same shape as y_con.
      integer, intent(in) :: y_cat(:, :) !! Categorical responses encoded 1..K, shape n by n_cat; masked entries may be zero.
      logical, intent(in) :: cat_observed(:, :) !! Observation mask for categorical responses, same shape as y_cat.
      integer, intent(in) :: n_levels(:) !! Number of levels K for each categorical response.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by q.
      integer, intent(in) :: n_iter !! Number of MCMC iterations; must be positive.
      real(dp), intent(in) :: prior_scale(:, :) !! p by p inverse-Wishart residual scale for all-continuous models.
      type(jomo1_result), intent(out) :: result !! Final imputation, parameter state, posterior means, and optional complete chains.
      real(dp), intent(in), optional :: beta_start(:, :) !! Optional q by p starting fixed-effect matrix; defaults to zero.
      real(dp), intent(in), optional :: omega_start(:, :) !! Optional p by p covariance start; default identity before constraints.
      logical, intent(in), optional :: store_chain !! If true, retain every beta and omega draw in result; defaults to false.
      integer :: n
      integer :: n_con
      integer :: n_cat
      integer :: p
      integer :: q
      integer :: iter
      integer :: i
      integer :: j
      integer :: info
      integer :: failures
      integer :: total_failures
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: omega(:, :)
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: xtx_inv(:, :)
      real(dp), allocatable :: bhat(:, :)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: scale(:, :)
      logical, allocatable :: observed(:, :)
      logical, allocatable :: fixed(:, :)
      logical :: keep_chain

      n = size(x, 1)
      n_con = size(y_con, 2)
      n_cat = size(y_cat, 2)
      if (size(y_con, 1) /= n .or. size(y_cat, 1) /= n) error stop "jomo1_mixed_mcmc: response/design row mismatch"
      if (size(n_levels) /= n_cat) error stop "jomo1_mixed_mcmc: n_levels size mismatch"
      if (n_iter <= 0) error stop "jomo1_mixed_mcmc: n_iter must be positive"
      p = latent_dimension(n_con, n_levels)
      q = size(x, 2)
      if (size(prior_scale, 1) /= p .or. size(prior_scale, 2) /= p) error stop "jomo1_mixed_mcmc: prior scale shape mismatch"
      if (n <= max(p, q)) error stop "jomo1_mixed_mcmc: need more rows than latent/design dimension"

      call initialize_latent_data(y_con, con_observed, y_cat, cat_observed, n_levels, result%latent, observed)
      allocate(beta(q, p), omega(p, p), fixed(p, p), xtx(q, q), xtx_inv(q, q), bhat(q, p))
      allocate(mean(n, p), residual(n, p), scale(p, p))
      beta = 0.0_dp
      omega = 0.0_dp
      do i = 1, p
         omega(i, i) = 1.0_dp
      end do
      if (present(beta_start)) then
         if (any(shape(beta_start) /= shape(beta))) error stop "jomo1_mixed_mcmc: beta_start shape mismatch"
         beta = beta_start
      end if
      if (present(omega_start)) then
         if (any(shape(omega_start) /= shape(omega))) error stop "jomo1_mixed_mcmc: omega_start shape mismatch"
         omega = omega_start
      end if
      call set_categorical_covariance(omega, n_con, n_levels, fixed)
      if (.not. is_spd(omega)) error stop "jomo1_mixed_mcmc: starting omega is not positive definite"

      xtx = matmul(transpose(x), x)
      call inverse_spd(xtx, xtx_inv, info)
      if (info /= 0) error stop "jomo1_mixed_mcmc: fixed-effect design is rank deficient"
      keep_chain = .false.
      if (present(store_chain)) keep_chain = store_chain
      if (keep_chain) then
         allocate(result%beta_chain(q, p, n_iter), result%omega_chain(p, p, n_iter))
      end if
      allocate(result%beta_mean(q, p), result%omega_mean(p, p))
      result%beta_mean = 0.0_dp
      result%omega_mean = 0.0_dp
      total_failures = 0

      do iter = 1, n_iter
         mean = matmul(x, beta)
         if (n_cat > 0) then
            call update_observed_categories(rng, y_cat, cat_observed, n_levels, n_con, mean, omega, result%latent, failures)
            total_failures = total_failures + failures
         end if

         bhat = matmul(xtx_inv, matmul(transpose(x), result%latent))
         call matrix_normal_sample(rng, bhat, xtx_inv, omega, beta, info)
         if (info /= 0) error stop "jomo1_mixed_mcmc: beta draw failed"
         mean = matmul(x, beta)
         residual = result%latent - mean

         if (n_cat == 0) then
            scale = prior_scale + matmul(transpose(residual), residual)
            call invwishart_sample(rng, real(n - 1, dp), scale, omega, info)
            if (info /= 0) error stop "jomo1_mixed_mcmc: inverse-Wishart residual draw failed"
         else
            call mh_covariance_update(rng, residual, fixed, omega)
         end if

         call impute_missing_rows(rng, observed, mean, omega, result%latent, info)
         if (info /= 0) error stop "jomo1_mixed_mcmc: conditional imputation failed"
         result%beta_mean = result%beta_mean + beta
         result%omega_mean = result%omega_mean + omega
         if (keep_chain) then
            result%beta_chain(:, :, iter) = beta
            result%omega_chain(:, :, iter) = omega
         end if
      end do

      result%beta_mean = result%beta_mean / real(n_iter, dp)
      result%omega_mean = result%omega_mean / real(n_iter, dp)
      allocate(result%beta(q, p), result%omega(p, p))
      result%beta = beta
      result%omega = omega
      allocate(result%continuous(n, n_con))
      if (n_con > 0) result%continuous = result%latent(:, 1:n_con)
      allocate(result%categorical(n, n_cat))
      if (n_cat > 0) call decode_categories(result%latent, n_con, n_levels, result%categorical)
      result%latent_rejection_failures = total_failures
      result%iterations = n_iter

      do j = 1, n_con
         do i = 1, n
            if (con_observed(i, j)) result%continuous(i, j) = y_con(i, j)
         end do
      end do
      do j = 1, n_cat
         do i = 1, n
            if (cat_observed(i, j)) result%categorical(i, j) = y_cat(i, j)
         end do
      end do
   end subroutine jomo1_mixed_mcmc

   subroutine jomo1con_mcmc(rng, y, observed, x, n_iter, prior_scale, result, beta_start, omega_start, store_chain)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all MCMC draws.
      real(dp), intent(in) :: y(:, :) !! Continuous response matrix, shape n by p; masked entries may contain any value.
      logical, intent(in) :: observed(:, :) !! Observation mask for y, true for observed values.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by q.
      integer, intent(in) :: n_iter !! Number of Gibbs iterations to run.
      real(dp), intent(in) :: prior_scale(:, :) !! Inverse-Wishart residual scale matrix, p by p.
      type(jomo1_result), intent(out) :: result !! Final continuous imputation, parameter state, means, and optional chains.
      real(dp), intent(in), optional :: beta_start(:, :) !! Optional q by p starting fixed-effect matrix.
      real(dp), intent(in), optional :: omega_start(:, :) !! Optional p by p starting residual covariance matrix.
      logical, intent(in), optional :: store_chain !! If true, retain every beta and omega draw.
      integer, allocatable :: y_cat(:, :)
      integer, allocatable :: n_levels(:)
      logical, allocatable :: cat_observed(:, :)

      allocate(y_cat(size(y, 1), 0), cat_observed(size(y, 1), 0), n_levels(0))
      call jomo1_mixed_mcmc(rng, y, observed, y_cat, cat_observed, n_levels, x, n_iter, prior_scale, result, &
         beta_start, omega_start, store_chain)
   end subroutine jomo1con_mcmc

   subroutine jomo1cat_mcmc(rng, y_cat, observed, n_levels, x, n_iter, prior_scale, result, beta_start, omega_start, store_chain)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all MCMC draws.
      integer, intent(in) :: y_cat(:, :) !! Categorical response matrix encoded 1..K, shape n by n_cat.
      logical, intent(in) :: observed(:, :) !! Observation mask for categorical responses.
      integer, intent(in) :: n_levels(:) !! Number of levels K for each categorical response.
      real(dp), intent(in) :: x(:, :) !! Fully observed fixed-effect design matrix, shape n by q.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale(:, :) !! p by p scale retained for API consistency; categorical covariance uses MH.
      type(jomo1_result), intent(out) :: result !! Final categories, latent values, parameter state, means, and optional chains.
      real(dp), intent(in), optional :: beta_start(:, :) !! Optional q by p starting fixed-effect matrix.
      real(dp), intent(in), optional :: omega_start(:, :) !! Optional p by p starting latent covariance matrix.
      logical, intent(in), optional :: store_chain !! If true, retain every beta and omega draw.
      real(dp), allocatable :: y_con(:, :)
      logical, allocatable :: con_observed(:, :)

      allocate(y_con(size(y_cat, 1), 0), con_observed(size(y_cat, 1), 0))
      call jomo1_mixed_mcmc(rng, y_con, con_observed, y_cat, observed, n_levels, x, n_iter, prior_scale, result, &
         beta_start, omega_start, store_chain)
   end subroutine jomo1cat_mcmc

   subroutine mh_covariance_update(rng, residual, fixed, omega)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for random-walk MH proposals and acceptance draws.
      real(dp), intent(in) :: residual(:, :) !! Current complete-data residual matrix, shape n by p.
      logical, intent(in) :: fixed(:, :) !! True for covariance elements fixed by latent-category identification constraints.
      real(dp), intent(inout) :: omega(:, :) !! Positive-definite covariance updated at each free upper-triangular element.
      integer :: n
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

      n = size(residual, 1)
      p = size(residual, 2)
      allocate(proposal(p, p))
      current_ll = gaussian_loglik(residual, omega)
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
            proposed_ll = gaussian_loglik(residual, proposal)
            if (log(rng_uniform(rng)) < proposed_ll - current_ll) then
               omega = proposal
               current_ll = proposed_ll
            end if
         end do
      end do
   end subroutine mh_covariance_update

   function gaussian_loglik(residual, omega) result(value)
      real(dp), intent(in) :: residual(:, :) !! Residual matrix with independent rows under covariance omega.
      real(dp), intent(in) :: omega(:, :) !! Symmetric positive-definite covariance matrix.
      real(dp) :: value
      real(dp), allocatable :: precision(:, :)
      integer :: i
      integer :: info

      allocate(precision(size(omega, 1), size(omega, 2)))
      call inverse_spd(omega, precision, info)
      if (info /= 0) then
         value = -huge(1.0_dp)
         return
      end if
      value = -0.5_dp * real(size(residual, 1), dp) * logdet_spd(omega)
      do i = 1, size(residual, 1)
         value = value - 0.5_dp * quadratic_form(residual(i, :), precision)
      end do
   end function gaussian_loglik

end module jomo_single_level
