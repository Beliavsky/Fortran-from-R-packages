! SPDX-License-Identifier: GPL-3.0-only
program test_pan_multivariate
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use pan, only : PAN_OK, dp, pan_mcmc, pan_prior, pan_result
   use pan_linalg, only : is_spd
   implicit none

   type(pan_prior) :: prior
   type(pan_result) :: fit
   integer :: i
   integer :: subj(15)
   integer :: xcol(2)
   integer :: zcol(2)
   real(dp) :: pred(15, 2)
   real(dp) :: y(15, 2)
   real(dp) :: original(15, 2)
   real(dp) :: nan

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   subj = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5]
   xcol = [1, 2]
   zcol = [1, 2]

   do i = 1, 15
      pred(i, 1) = 1.0_dp
      pred(i, 2) = -1.0_dp + real(mod(i - 1, 3), dp)
      y(i, 1) = 0.5_dp + 0.4_dp * pred(i, 2) + 0.08_dp * real(i, dp)
      y(i, 2) = -0.2_dp + 0.7_dp * pred(i, 2) + 0.05_dp * real(i, dp)
   end do
   original = y
   y(2, 1) = nan
   y(7, 2) = nan
   y(12, :) = nan

   prior%a = 2.0_dp
   allocate(prior%binv(2, 2))
   prior%binv = reshape([1.0_dp, 0.2_dp, 0.2_dp, 1.2_dp], [2, 2])
   prior%c = 4.0_dp
   allocate(prior%dinv(4, 4))
   prior%dinv = 0.0_dp
   do i = 1, 4
      prior%dinv(i, i) = 2.0_dp
   end do

   call pan_mcmc(y, subj, pred, xcol, zcol, prior, 777, 6, fit)

   call assert_true(fit%status == PAN_OK, "multivariate pan sampler status")
   call assert_true(all(shape(fit%beta) == [2, 2, 6]), "multivariate beta chain shape")
   call assert_true(all(shape(fit%psi) == [4, 4, 6]), "full random covariance chain shape")
   call assert_true(is_spd(fit%last%sigma), "multivariate residual covariance SPD")
   call assert_true(is_spd(fit%last%psi), "multivariate random covariance SPD")
   call assert_true(.not. any(ieee_is_nan(fit%y)), "all response entries completed")
   call assert_close(fit%y(1, 1), original(1, 1), 0.0_dp, "observed component one preserved")
   call assert_close(fit%y(7, 1), original(7, 1), 0.0_dp, "paired observed component preserved")
   call assert_close(fit%y(2, 2), original(2, 2), 0.0_dp, "other paired observed component preserved")

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

end program test_pan_multivariate
