! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_inference
   use sde
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   integer, parameter :: n_ar = 401
   integer, parameter :: n_ou = 161
   real(dp), parameter :: dt = 0.05_dp
   real(dp), parameter :: true_ou(3) = [0.9_dp, 0.7_dp, 0.35_dp]
   real(dp) :: ar_series(n_ar), equations1(1), generator_equations(1)
   real(dp), allocatable :: ou_matrix(:, :), ou_series(:)
   type(estimating_result) :: simple_fit, generator_fit, martingale_fit
   type(gmm_result) :: gmm_fit
   type(sde_aic_result) :: aic_fit
   type(sde_divergence_result) :: divergence, fitted_divergence
   real(dp) :: exact_ll, euler_ll, dc_ll, criterion, target_intercept
   integer :: i, failures

   failures = 0
   call seed_rng(271828182_i64)

   ar_series(1) = 0.0_dp
   do i = 2, n_ar
      ar_series(i) = 0.65_dp*ar_series(i-1)+random_normal(sd=0.4_dp)
   end do

   call evaluate_simple_estimating(ar_series, [0.65_dp], 1, ar_estimating, equations1)
   call assert_true(abs(equations1(1)) < 20.0_dp, "simple estimating equation evaluation", failures)
   call fit_simple_estimating(ar_series, 1, ar_estimating, [0.2_dp], simple_fit, &
      lower=[-0.99_dp], upper=[0.99_dp], max_iterations=600)
   call assert_true(abs(simple_fit%estimate(1)-0.65_dp) < 0.12_dp, &
      "simple estimating function fit", failures)
   call assert_true(simple_fit%objective < 1.0e-5_dp, &
      "simple estimating objective", failures)

   call fit_linear_martingale(ar_series, ar_conditional_mean, ar_conditional_variance, &
      ar_weights, [0.2_dp], martingale_fit, lower=[-0.99_dp], upper=[0.99_dp], &
      max_iterations=600)
   call assert_true(abs(martingale_fit%estimate(1)-0.65_dp) < 0.12_dp, &
      "linear martingale fit", failures)
   call assert_true(martingale_fit%objective < 1.0e-5_dp, &
      "linear martingale objective", failures)

   call fit_gmm(ar_series, 1.0_dp, 2, ar_moments, [0.2_dp, 0.3_dp], gmm_fit, &
      lower=[-0.99_dp, 0.02_dp], upper=[0.99_dp, 2.0_dp], max_stage2_iterations=5, &
      parameter_tolerance=2.0e-3_dp, objective_tolerance=1.0e-7_dp, max_lag=4, &
      optimizer_iterations=700)
   call assert_true(abs(gmm_fit%estimate(1)-0.65_dp) < 0.15_dp, "GMM AR coefficient", failures)
   call assert_true(abs(gmm_fit%estimate(2)-0.16_dp) < 0.10_dp, "GMM innovation variance", failures)
   call assert_true(all(ieee_is_finite(gmm_fit%hessian)), "GMM Hessian finite", failures)

   call simulate_ou_exact([0.2_dp], dt, n_ou-1, true_ou, ou_matrix)
   allocate(ou_series(n_ou))
   ou_series = ou_matrix(:, 1)

   exact_ll = ou_log_likelihood(ou_series, dt, true_ou)
   euler_ll = euler_log_likelihood(ou_series, dt, true_ou, ou_coefficient_drift, &
      ou_coefficient_diffusion)
   dc_ll = dc_log_likelihood(ou_series, dt, true_ou, ou_state_drift, ou_state_diffusion, &
      ou_state_drift_x, zero_state, zero_state)
   criterion = sde_aic(ou_series, dt, true_ou, ou_state_drift, ou_state_diffusion, &
      ou_state_drift_x, zero_state, zero_state)
   call assert_true(ieee_is_finite(exact_ll), "exact OU log likelihood finite", failures)
   call assert_true(ieee_is_finite(euler_ll), "Euler log likelihood finite", failures)
   call assert_true(ieee_is_finite(dc_ll), "Dacunha-Castelle log likelihood finite", failures)
   call assert_close(criterion, -2.0_dp*dc_ll+6.0_dp, 2.0e-11_dp, "SDE AIC formula", failures)

   call fit_sde_aic(ou_series, dt, intercept_drift_state, fixed_diffusion_state, &
      fixed_drift_x_state, zero_state, zero_state, [0.4_dp], aic_fit, lower=[-1.0_dp], &
      upper=[2.0_dp], max_iterations=500)
   call assert_true(size(aic_fit%estimate) == 1, "SDE AIC fit dimensions", failures)
   call assert_true(ieee_is_finite(aic_fit%aic), "SDE AIC fit finite", failures)

   target_intercept = 0.7_dp*sum(ou_series)/real(size(ou_series), dp)
   call evaluate_generator_estimating(ou_series, [target_intercept], 1, &
      fixed_intercept_coefficient_drift, fixed_coefficient_diffusion, linear_h_derivatives, &
      generator_equations)
   call assert_true(abs(generator_equations(1)) < 1.0e-11_dp, &
      "generator equation root", failures)
   call fit_generator_estimating(ou_series, 1, fixed_intercept_coefficient_drift, &
      fixed_coefficient_diffusion, linear_h_derivatives, [0.2_dp], generator_fit, &
      lower=[-1.0_dp], upper=[2.0_dp], max_iterations=500)
   call assert_true(abs(generator_fit%estimate(1)-target_intercept) < 2.0e-5_dp, &
      "generator estimating fit", failures)

   call sde_divergence_test(ou_series, dt, true_ou, [0.75_dp, 0.7_dp, 0.35_dp], &
      ou_state_drift, ou_state_diffusion, ou_state_drift_x, zero_state, zero_state, &
      divergence, n_simulations=2500)
   call assert_true(divergence%likelihood_ratio >= 0.0_dp, "divergence LRT nonnegative", failures)
   call assert_probability(divergence%p_likelihood_ratio, "divergence LRT p-value", failures)
   call assert_probability(divergence%p_divergence, "divergence p-value", failures)

   call fit_and_test_sde_divergence(ou_series, dt, [0.75_dp], [0.4_dp], &
      intercept_drift_state, fixed_diffusion_state, fixed_drift_x_state, zero_state, zero_state, &
      fitted_divergence, n_simulations=1000, lower=[-1.0_dp], upper=[2.0_dp], max_iterations=500)
   call assert_true(size(fitted_divergence%theta1) == 1, &
      "fitted divergence parameter dimensions", failures)
   call assert_probability(fitted_divergence%p_likelihood_ratio, &
      "fitted divergence LRT p-value", failures)

   if (failures /= 0) then
      write(*, '(a, i0)') "test_inference: failures = ", failures
      error stop 1
   end if
   write(*, '(a)') "test_inference: PASS"

contains

   subroutine ar_estimating(y, x, theta, values)
      real(dp), intent(in) :: y, x, theta(:)
      real(dp), intent(out) :: values(:)
      values(1) = x*(y-theta(1)*x)
   end subroutine ar_estimating

   pure function ar_conditional_mean(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = theta(1)*x
   end function ar_conditional_mean

   pure function ar_conditional_variance(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = 0.16_dp+0.0_dp*(x+sum(theta))
   end function ar_conditional_variance

   subroutine ar_weights(order, index, x, theta, weight)
      integer, intent(in) :: order, index
      real(dp), intent(in) :: x, theta(:)
      real(dp), intent(out) :: weight
      if (order == 1 .and. index == 1) then
         weight = x
      else
         weight = 0.0_dp
      end if
      weight = weight+0.0_dp*sum(theta)
   end subroutine ar_weights

   subroutine ar_moments(y, x, theta, local_dt, moments)
      real(dp), intent(in) :: y, x, theta(:), local_dt
      real(dp), intent(out) :: moments(:)
      real(dp) :: residual
      residual = y-theta(1)*x
      moments(1) = x*residual
      moments(2) = residual*residual-theta(2)
      moments = moments+0.0_dp*local_dt
   end subroutine ar_moments

   pure function ou_coefficient_drift(t, x, theta) result(value)
      real(dp), intent(in) :: t, x, theta(:)
      real(dp) :: value
      value = theta(1)-theta(2)*x+0.0_dp*t
   end function ou_coefficient_drift

   pure function ou_coefficient_diffusion(t, x, theta) result(value)
      real(dp), intent(in) :: t, x, theta(:)
      real(dp) :: value
      value = theta(3)+0.0_dp*(t+x)
   end function ou_coefficient_diffusion

   pure function ou_state_drift(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = theta(1)-theta(2)*x
   end function ou_state_drift

   pure function ou_state_diffusion(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = theta(3)+0.0_dp*x
   end function ou_state_diffusion

   pure function ou_state_drift_x(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = -theta(2)+0.0_dp*x
   end function ou_state_drift_x

   pure function intercept_drift_state(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = theta(1)-0.7_dp*x
   end function intercept_drift_state

   pure function fixed_diffusion_state(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = 0.35_dp+0.0_dp*(x+sum(theta))
   end function fixed_diffusion_state

   pure function fixed_drift_x_state(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = -0.7_dp+0.0_dp*(x+sum(theta))
   end function fixed_drift_x_state

   pure function zero_state(x, theta) result(value)
      real(dp), intent(in) :: x, theta(:)
      real(dp) :: value
      value = 0.0_dp*(x+real(size(theta), dp))
   end function zero_state

   pure function fixed_intercept_coefficient_drift(t, x, theta) result(value)
      real(dp), intent(in) :: t, x, theta(:)
      real(dp) :: value
      value = theta(1)-0.7_dp*x+0.0_dp*t
   end function fixed_intercept_coefficient_drift

   pure function fixed_coefficient_diffusion(t, x, theta) result(value)
      real(dp), intent(in) :: t, x, theta(:)
      real(dp) :: value
      value = 0.35_dp+0.0_dp*(t+x+sum(theta))
   end function fixed_coefficient_diffusion

   subroutine linear_h_derivatives(index, x, theta, h_x, h_xx)
      integer, intent(in) :: index
      real(dp), intent(in) :: x, theta(:)
      real(dp), intent(out) :: h_x, h_xx
      h_x = 1.0_dp
      h_xx = 0.0_dp
      if (index /= 1) then
         h_x = 0.0_dp
         h_xx = 0.0_dp
      end if
      h_x = h_x+0.0_dp*(x+sum(theta))
   end subroutine linear_h_derivatives

   subroutine assert_close(actual, target, tolerance, label, count)
      real(dp), intent(in) :: actual, target, tolerance
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. ieee_is_finite(actual) .or. abs(actual-target) > tolerance) then
         write(*, '(a, 2(1x, es23.15), a, es12.4)') "FAIL "//trim(label)//":", actual, target, &
            " tolerance=", tolerance
         count = count+1
      end if
   end subroutine assert_close

   subroutine assert_probability(value, label, count)
      real(dp), intent(in) :: value
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      call assert_true(ieee_is_finite(value) .and. value >= 0.0_dp .and. value <= 1.0_dp, &
         label, count)
   end subroutine assert_probability

   subroutine assert_true(condition, label, count)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. condition) then
         write(*, '(a)') "FAIL "//trim(label)
         count = count+1
      end if
   end subroutine assert_true

end program test_inference
