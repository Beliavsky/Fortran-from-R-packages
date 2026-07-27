! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_models
   use sde
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   real(dp), parameter :: ou_theta(3) = [1.2_dp, 0.8_dp, 0.5_dp]
   real(dp), parameter :: gbm_theta(2) = [0.07_dp, 0.25_dp]
   real(dp), parameter :: cir_theta(3) = [0.9_dp, 1.4_dp, 0.45_dp]
   real(dp), parameter :: probabilities(5) = [0.001_dp, 0.1_dp, 0.5_dp, 0.9_dp, 0.999_dp]
   real(dp) :: p, q, x, expected, sample
   integer :: i, failures

   failures = 0

   do i = 1, size(probabilities)
      p = probabilities(i)
      q = normal_quantile(p)
      call assert_close(normal_cdf(q), p, 2.0e-12_dp, "normal quantile round trip", failures)

      q = gamma_quantile(p, 2.5_dp, 1.7_dp)
      call assert_close(gamma_cdf(q, 2.5_dp, 1.7_dp), p, 3.0e-11_dp, &
         "gamma quantile round trip", failures)

      q = noncentral_chi_square_quantile(p, 5.0_dp, 2.0_dp)
      call assert_close(noncentral_chi_square_cdf(q, 5.0_dp, 2.0_dp), p, 3.0e-10_dp, &
         "noncentral chi-square quantile round trip", failures)
   end do

   expected = ou_theta(1)/ou_theta(2)+(0.3_dp-ou_theta(1)/ou_theta(2))*exp(-ou_theta(2)*0.4_dp)
   call assert_close(ou_conditional_mean(0.4_dp, 0.3_dp, ou_theta), expected, 2.0e-14_dp, &
      "OU conditional mean", failures)
   expected = ou_theta(3)**2*(1.0_dp-exp(-2.0_dp*ou_theta(2)*0.4_dp))/(2.0_dp*ou_theta(2))
   call assert_close(ou_conditional_variance(0.4_dp, ou_theta), expected, 2.0e-14_dp, &
      "OU conditional variance", failures)

   q = ou_conditional_quantile(0.73_dp, 0.4_dp, 0.3_dp, ou_theta)
   call assert_close(ou_conditional_cdf(q, 0.4_dp, 0.3_dp, ou_theta), 0.73_dp, 2.0e-12_dp, &
      "OU conditional distribution", failures)
   q = ou_stationary_quantile(0.37_dp, ou_theta)
   call assert_close(ou_stationary_cdf(q, ou_theta), 0.37_dp, 2.0e-12_dp, &
      "OU stationary distribution", failures)

   q = gbm_conditional_quantile(0.81_dp, 0.75_dp, 10.0_dp, gbm_theta)
   call assert_close(gbm_conditional_cdf(q, 0.75_dp, 10.0_dp, gbm_theta), 0.81_dp, 2.0e-12_dp, &
      "GBM conditional distribution", failures)
   call assert_true(gbm_conditional_pdf(q, 0.75_dp, 10.0_dp, gbm_theta) > 0.0_dp, &
      "GBM density is positive", failures)

   q = cir_conditional_quantile(0.64_dp, 0.3_dp, 0.7_dp, cir_theta)
   call assert_close(cir_conditional_cdf(q, 0.3_dp, 0.7_dp, cir_theta), 0.64_dp, 2.0e-9_dp, &
      "CIR conditional distribution", failures)
   q = cir_stationary_quantile(0.42_dp, cir_theta)
   call assert_close(cir_stationary_cdf(q, cir_theta), 0.42_dp, 2.0e-11_dp, &
      "CIR stationary distribution", failures)

   call seed_rng(8675309_i64)
   sample = ou_conditional_random(0.1_dp, 0.2_dp, ou_theta)
   call assert_true(ieee_is_finite(sample), "OU random draw is finite", failures)
   sample = gbm_conditional_random(0.1_dp, 1.0_dp, gbm_theta)
   call assert_true(sample > 0.0_dp, "GBM random draw is positive", failures)
   sample = cir_conditional_random(0.1_dp, 0.4_dp, cir_theta)
   call assert_true(sample >= 0.0_dp, "CIR random draw is nonnegative", failures)
   sample = cir_stationary_random(cir_theta)
   call assert_true(sample >= 0.0_dp, "stationary CIR random draw is nonnegative", failures)

   x = integrate_adaptive(square_function, 0.0_dp, 1.0_dp, [real(dp) ::])
   call assert_close(x, 1.0_dp/3.0_dp, 1.0e-10_dp, "adaptive integration", failures)

   if (failures /= 0) then
      write(*, '(a, i0)') "test_models: failures = ", failures
      error stop 1
   end if
   write(*, '(a)') "test_models: PASS"

contains

   pure function square_function(value, theta) result(answer)
      real(dp), intent(in) :: value, theta(:)
      real(dp) :: answer
      answer = value*value+0.0_dp*real(size(theta), dp)
   end function square_function

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

   subroutine assert_true(condition, label, count)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. condition) then
         write(*, '(a)') "FAIL "//trim(label)
         count = count+1
      end if
   end subroutine assert_true

end program test_models
