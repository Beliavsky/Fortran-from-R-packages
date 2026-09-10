! SPDX-License-Identifier: GPL-2.0-or-later
program test_stats_threshold
   use wavethresh
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: x(4)
   real(dp) :: sigma(4)
   real(dp) :: input(8)
   real(dp), allocatable :: fourier(:)
   real(dp), allocatable :: inverse(:)
   type(wd_t) :: object
   type(wd_t) :: shrunk
   integer :: levels(2)

   x = [0.5_dp, -1.0_dp, 2.0_dp, -3.0_dp]
   sigma = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
   call assert_close(sure(x), 0.5_dp, tol, "sure reference")
   call assert_close(newsure(sigma, x), 0.5_dp, tol, "newsure reference")
   call assert_close(shannon_entropy([0.5_dp, 0.5_dp]), 0.6931471805599453_dp, tol, "entropy reference")
   call assert_close(ssq([1.0_dp, 2.0_dp], [2.0_dp, 4.0_dp]), 5.0_dp, tol, "ssq reference")
   call assert_close(dclaw(0.0_dp), 0.5984163940411785_dp, 2.0e-12_dp, "claw density reference")
   call assert_close(pclaw(0.0_dp), 0.5_dp, 2.0e-12_dp, "claw CDF symmetry")

   input = [0.0_dp, 1.0_dp, 4.0_dp, -2.0_dp, 3.0_dp, 0.5_dp, -1.0_dp, 2.0_dp]
   fourier = rfft(input)
   inverse = rfftinv(fourier)
   call assert_close(maxval(abs(inverse - input)), 0.0_dp, 3.0e-12_dp, "rfft round trip")

   object = wd(input, filter_number=1.0_dp, family="DaubExPhase")
   levels = [1, 2]
   shrunk = threshold_wd(object, levels=levels, threshold_type="hard", policy="manual", value=huge(1.0_dp) / 4.0_dp)
   call assert_close(maxval(abs(shrunk%detail(1)%values)), 0.0_dp, tol, "manual threshold level 1")
   call assert_close(maxval(abs(shrunk%detail(2)%values)), 0.0_dp, tol, "manual threshold level 2")
   call assert_true(dof(shrunk) <= dof(object), "threshold reduces degrees of freedom")

   print *, "test_stats_threshold: PASS"
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must hold.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Value produced by the translation.
      real(dp), intent(in) :: expected !! Deterministic reference value.
      real(dp), intent(in) :: tolerance !! Maximum accepted absolute error.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close
end program test_stats_threshold
