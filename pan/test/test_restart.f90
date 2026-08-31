! SPDX-License-Identifier: GPL-3.0-only
program test_restart
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use pan, only : PAN_OK, dp, pan_mcmc, pan_prior, pan_result
   implicit none

   type(pan_prior) :: prior
   type(pan_result) :: first_run
   type(pan_result) :: restart_a
   type(pan_result) :: restart_b
   integer :: i
   integer :: subj(9)
   integer :: xcol(2)
   integer :: zcol(1)
   real(dp) :: pred(9, 2)
   real(dp) :: y(9, 1)
   real(dp) :: nan

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   subj = [1, 1, 1, 2, 2, 2, 3, 3, 3]
   xcol = [1, 2]
   zcol = [1]

   do i = 1, 9
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(mod(i - 1, 3), dp)
      y(i, 1) = 0.7_dp + 0.3_dp * pred(i, 2) + 0.04_dp * real(i, dp)
   end do
   y(5, 1) = nan

   prior%a = 1.0_dp
   prior%binv = reshape([1.0_dp], [1, 1])
   prior%c = 1.0_dp
   prior%dinv = reshape([1.0_dp], [1, 1])

   call pan_mcmc(y, subj, pred, xcol, zcol, prior, 11, 3, first_run)
   call assert_true(first_run%status == PAN_OK, "initial chain for restart succeeds")

   call pan_mcmc(y, subj, pred, xcol, zcol, prior, 919, 2, restart_a, start=first_run%last)
   call pan_mcmc(y, subj, pred, xcol, zcol, prior, 919, 2, restart_b, start=first_run%last)

   call assert_true(restart_a%status == PAN_OK, "restart chain succeeds")
   call assert_true(restart_b%status == PAN_OK, "replicated restart chain succeeds")
   call assert_true(.not. ieee_is_nan(restart_a%y(5, 1)), "restart chain keeps missing value imputed")
   call assert_close(restart_a%beta(1, 1, 1), restart_b%beta(1, 1, 1), 0.0_dp, "restart is deterministic for a seed")
   call assert_close(restart_a%sigma(1, 1, 2), restart_b%sigma(1, 1, 2), 0.0_dp, "restart covariance reproducible")
   call assert_close(restart_a%y(1, 1), y(1, 1), 0.0_dp, "restart preserves observed data")

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

end program test_restart
