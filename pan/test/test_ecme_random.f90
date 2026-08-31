! SPDX-License-Identifier: GPL-3.0-only
program test_ecme_random
   use pan, only : PAN_OK, dp, ecme_fit, ecme_result
   use pan_linalg, only : is_spd
   implicit none

   type(ecme_result) :: fit
   integer :: i
   integer :: occ(30)
   integer :: subj(30)
   integer :: xcol(2)
   integer :: zcol(1)
   real(dp) :: pred(30, 2)
   real(dp) :: y(30)
   real(dp), parameter :: btrue(10) = [ &
      -0.8_dp, -0.4_dp, -0.2_dp, 0.0_dp, 0.15_dp, 0.3_dp, 0.45_dp, 0.6_dp, -0.1_dp, 0.25_dp ]

   xcol = [1, 2]
   zcol = [1]

   do i = 1, 30
      subj(i) = (i - 1) / 3 + 1
      occ(i) = mod(i - 1, 3) + 1
      pred(i, 1) = 1.0_dp
      pred(i, 2) = real(occ(i) - 1, dp)
      y(i) = 2.0_dp + 0.5_dp * pred(i, 2) + btrue(subj(i)) + 0.06_dp * real(mod(i, 3) - 1, dp)
   end do

   call ecme_fit(y, subj, occ, pred, xcol, fit, zcol=zcol, maxits=2000, tol=1.0e-7_dp)

   call assert_true(fit%status == PAN_OK, "random-intercept ecme status")
   call assert_true(fit%converged, "random-intercept ecme converged")
   call assert_true(fit%sigma2 > 0.0_dp, "random-intercept residual variance positive")
   call assert_true(is_spd(fit%psi), "random-intercept psi positive definite")
   call assert_true(all(shape(fit%bhat) == [1, 10]), "empirical Bayes random-effect shape")
   call assert_true(all(shape(fit%cov_b) == [1, 1, 10]), "random-effect covariance shape")
   call assert_true(fit%iter >= 2, "iterative path used")
   call assert_true(fit%loglik(fit%iter) >= fit%loglik(1) - 1.0e-8_dp, "marginal likelihood does not decrease overall")

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

end program test_ecme_random
