! SPDX-License-Identifier: GPL-3.0-only
program test_errors
   use pan, only : PAN_ERR_ARGUMENT, dp, pan_mcmc, pan_prior, pan_result
   implicit none

   type(pan_prior) :: prior
   type(pan_result) :: fit
   integer :: subj(4)
   integer :: xcol(1)
   integer :: zcol(1)
   real(dp) :: pred(4, 1)
   real(dp) :: y(4, 1)

   subj = [1, 2, 1, 2]
   pred = 1.0_dp
   y(:, 1) = [1.0_dp, 2.0_dp, 1.5_dp, 2.5_dp]
   xcol = [1]
   zcol = [1]
   prior%a = 1.0_dp
   prior%binv = reshape([1.0_dp], [1, 1])
   prior%c = 1.0_dp
   prior%dinv = reshape([1.0_dp], [1, 1])

   call pan_mcmc(y, subj, pred, xcol, zcol, prior, 1, 2, fit)
   call assert_true(fit%status == PAN_ERR_ARGUMENT, "unsorted subject labels rejected")

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


   subroutine finish_tests()
   end subroutine finish_tests

end program test_errors
