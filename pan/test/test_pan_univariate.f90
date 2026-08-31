! SPDX-License-Identifier: GPL-3.0-only
program test_pan_univariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use pan, only : PAN_OK, dp, pan_mcmc, pan_prior, pan_result
   use pan_linalg, only : is_spd
   implicit none

   type(pan_prior) :: prior
   type(pan_result) :: fit
   integer :: i
   integer :: subj(12)
   integer :: xcol(2)
   integer :: zcol(1)
   real(dp) :: pred(12, 2)
   real(dp) :: y(12, 1)
   real(dp) :: nan
   real(dp), parameter :: observed_y(12) = [ &
      1.1_dp, 1.4_dp, 1.7_dp, 1.3_dp, 1.6_dp, 1.9_dp, &
      0.8_dp, 1.1_dp, 1.5_dp, 1.5_dp, 1.8_dp, 2.1_dp ]

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   subj = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
   xcol = [1, 2]
   zcol = [1]

   do i = 1, 12
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(mod(i - 1, 3), dp)
      y(i, 1) = observed_y(i)
   end do
   y(2, 1) = nan
   y(8, 1) = nan

   prior%a = 1.0_dp
   prior%binv = reshape([1.0_dp], [1, 1])
   prior%c = 1.0_dp
   prior%dinv = reshape([1.0_dp], [1, 1])

   call pan_mcmc(y, subj, pred, xcol, zcol, prior, 2468, 8, fit)

   call assert_true(fit%status == PAN_OK, "univariate pan sampler status")
   call assert_true(all(shape(fit%beta) == [2, 1, 8]), "beta chain shape")
   call assert_true(all(shape(fit%sigma) == [1, 1, 8]), "sigma chain shape")
   call assert_true(all(shape(fit%psi) == [1, 1, 8]), "psi chain shape")
   call assert_true(fit%sigma(1, 1, 8) > 0.0_dp, "final residual variance positive")
   call assert_true(fit%psi(1, 1, 8) > 0.0_dp, "final random-effect variance positive")
   call assert_true(.not. ieee_is_nan(fit%y(2, 1)), "first missing value imputed")
   call assert_true(.not. ieee_is_nan(fit%y(8, 1)), "second missing value imputed")

   do i = 1, 12
      if (i /= 2 .and. i /= 8) call assert_close(fit%y(i, 1), observed_y(i), 0.0_dp, "observed y preserved")
   end do

   call assert_true(is_spd(fit%last%sigma), "restart residual covariance is SPD")
   call assert_true(is_spd(fit%last%psi), "restart random-effect covariance is SPD")
   call finish_tests()

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Test predicate expected to be true.
      character(len=*), intent(in) :: message !! Diagnostic printed if the predicate is false.

      if (.not. condition) then
         write (*, '(a)') "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, message)
      real(dp), intent(in) :: actual !! Computed scalar value.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tol !! Maximum allowed absolute error.
      character(len=*), intent(in) :: message !! Diagnostic printed if the comparison fails.

      if (abs(actual - expected) > tol) then
         write (*, '(a,2(es16.8,1x))') "FAIL: " // trim(message) // " actual/ref: ", actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine finish_tests()
   end subroutine finish_tests

end program test_pan_univariate
