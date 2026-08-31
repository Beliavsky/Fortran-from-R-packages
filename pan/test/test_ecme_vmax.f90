! SPDX-License-Identifier: GPL-3.0-only
program test_ecme_vmax
   use pan, only : PAN_OK, dp, ecme_fit, ecme_result
   implicit none

   type(ecme_result) :: fit
   integer :: i
   integer :: occ(6)
   integer :: subj(6)
   integer :: xcol(2)
   real(dp) :: pred(6, 2)
   real(dp) :: vmax(2, 2)
   real(dp) :: y(6)

   subj = [1, 1, 2, 2, 3, 3]
   occ = [1, 2, 1, 2, 1, 2]
   xcol = [1, 2]
   y = [1.0_dp, 1.9_dp, 1.4_dp, 2.7_dp, 1.7_dp, 3.6_dp]

   do i = 1, 6
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(i - 1, dp)
   end do

   vmax = reshape([1.0_dp, 0.4_dp, 0.4_dp, 1.0_dp], [2, 2])
   call ecme_fit(y, subj, occ, pred, xcol, fit, vmax=vmax)

   call assert_true(fit%status == PAN_OK, "correlated-vmax GLS status")
   call assert_close(fit%beta(1), 0.8213675213675214_dp, 1.0e-12_dp, "correlated-vmax intercept")
   call assert_close(fit%beta(2), 0.4914529914529915_dp, 1.0e-12_dp, "correlated-vmax slope")
   call assert_close(fit%sigma2, 0.4593542260208928_dp, 1.0e-12_dp, "correlated-vmax residual scale")

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

end program test_ecme_vmax
