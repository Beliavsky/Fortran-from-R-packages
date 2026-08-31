program test_substantive_models
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_substantive, only : normal_cdf, linear_loglik, binary_probit_loglik, ordinal_probit_loglik
   use jomo_substantive, only : cox_partial_loglik_ordered, sample_binary_probit_latent
   use jomo_substantive, only : sample_ordinal_probit_latent, update_ordinal_thresholds
   use jomo_substantive, only : sample_gaussian_coefficients, sample_linear_variance, cox_coordinate_newton
   implicit none

   integer, parameter :: n = 12
   real(dp) :: x(n, 2)
   real(dp) :: y_linear(n)
   real(dp) :: beta(2)
   real(dp) :: beta_draw(2)
   real(dp) :: latent(n)
   real(dp) :: residual(n)
   real(dp) :: variance
   integer :: y_binary(n)
   integer :: y_ordinal(n)
   real(dp) :: thresholds(2)
   logical :: event(n)
   logical :: converged
   type(rng_state) :: rng
   integer :: i
   integer :: info
   real(dp) :: ll0
   real(dp) :: ll1

   do i = 1, n
      x(i, 1) = 1.0_dp
      x(i, 2) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
   end do
   beta = [0.2_dp, 0.9_dp]
   y_linear = matmul(x, beta) + 0.05_dp * [(sin(real(i, dp)), i = 1, n)]
   if (.not. ieee_is_finite(linear_loglik(y_linear, x, beta, 0.5_dp))) error stop "linear likelihood failed"
   if (abs(normal_cdf(0.0_dp) - 0.5_dp) > 1.0e-14_dp) error stop "normal CDF failed"

   do i = 1, n
      if (dot_product(x(i, :), beta) < 0.0_dp) then
         y_binary(i) = 1
      else
         y_binary(i) = 2
      end if
   end do
   if (.not. ieee_is_finite(binary_probit_loglik(y_binary, x, beta))) error stop "binary likelihood failed"
   call rng_seed(rng, 1234567_i8)
   call sample_binary_probit_latent(rng, y_binary, x, beta, latent)
   if (any((y_binary == 1) .and. (latent >= 0.0_dp))) error stop "binary latent sign failed"
   if (any((y_binary == 2) .and. (latent <= 0.0_dp))) error stop "binary latent sign failed"
   call sample_gaussian_coefficients(rng, latent, x, 1.0_dp, beta_draw, info)
   if (info /= 0 .or. .not. all(ieee_is_finite(beta_draw))) error stop "Gaussian coefficient draw failed"

   thresholds = [-0.4_dp, 0.7_dp]
   do i = 1, n
      if (dot_product(x(i, :), beta) < thresholds(1)) then
         y_ordinal(i) = 1
      else if (dot_product(x(i, :), beta) < thresholds(2)) then
         y_ordinal(i) = 2
      else
         y_ordinal(i) = 3
      end if
   end do
   if (.not. ieee_is_finite(ordinal_probit_loglik(y_ordinal, x, beta, thresholds))) &
      error stop "ordinal likelihood failed"
   call sample_ordinal_probit_latent(rng, y_ordinal, x, beta, thresholds, latent)
   do i = 1, n
      select case (y_ordinal(i))
      case (1)
         if (latent(i) >= thresholds(1)) error stop "ordinal latent interval failed"
      case (2)
         if (latent(i) <= thresholds(1) .or. latent(i) >= thresholds(2)) error stop "ordinal latent interval failed"
      case (3)
         if (latent(i) <= thresholds(2)) error stop "ordinal latent interval failed"
      end select
   end do
   call update_ordinal_thresholds(rng, y_ordinal, latent, thresholds)
   if (thresholds(1) >= thresholds(2)) error stop "ordinal threshold ordering failed"

   event = .false.
   event(2) = .true.
   event(4) = .true.
   event(7) = .true.
   event(10) = .true.
   beta = [0.0_dp, 0.0_dp]
   ll0 = cox_partial_loglik_ordered(event, x, beta)
   call cox_coordinate_newton(event, x, beta, converged=converged)
   ll1 = cox_partial_loglik_ordered(event, x, beta)
   if (.not. ieee_is_finite(ll1) .or. ll1 < ll0 - 1.0e-10_dp) error stop "Cox Newton step reduced likelihood"

   residual = y_linear - matmul(x, [0.2_dp, 0.9_dp])
   call sample_linear_variance(rng, residual, 1.0_dp, variance)
   if (.not. ieee_is_finite(variance) .or. variance <= 0.0_dp) error stop "linear variance draw failed"

   print '(a)', "test_substantive_models: PASS"
end program test_substantive_models
