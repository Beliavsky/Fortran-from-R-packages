! SPDX-License-Identifier: GPL-3.0-only
program test_ecme_gls
   use pan, only : PAN_OK, dp, ecme_fit, ecme_result
   implicit none

   type(ecme_result) :: fit
   integer :: i
   integer :: occ(8)
   integer :: subj(8)
   integer :: xcol(2)
   real(dp) :: pred(8, 2)
   real(dp) :: y(8)

   subj = [1, 1, 2, 2, 3, 3, 4, 4]
   occ = [1, 2, 1, 2, 1, 2, 1, 2]
   xcol = [1, 2]

   do i = 1, 8
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(i - 1, dp)
      y(i) = 1.5_dp + 0.75_dp * pred(i, 2) + 0.1_dp * (-1.0_dp)**i
   end do

   call ecme_fit(y, subj, occ, pred, xcol, fit)

   call assert_true(fit%status == PAN_OK, "GLS ecme status")
   call assert_true(fit%converged, "GLS path converged")
   call assert_close(fit%beta(1), 1.4666666666666668_dp, 1.0e-12_dp, "GLS intercept")
   call assert_close(fit%beta(2), 0.7595238095238095_dp, 1.0e-12_dp, "GLS slope")
   call assert_close(fit%sigma2, 0.009523809523809518_dp, 1.0e-12_dp, "GLS ML residual variance")
   call assert_true(size(fit%psi, 1) == 0, "GLS path has no random-effect covariance")
   call assert_true(fit%cov_beta(1, 1) > 0.0_dp, "GLS beta covariance positive")

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

end program test_ecme_gls
