! SPDX-License-Identifier: GPL-2.0-or-later
program test_bayes
   use ieee_arithmetic, only : ieee_is_finite
   use wavethresh_types, only : dp
   use wavethresh_bayes, only : bayes_thr
   implicit none

   real(dp) :: data(32)
   real(dp) :: zeros(32)
   real(dp), allocatable :: estimate(:)
   real(dp), allocatable :: zero_estimate(:)
   integer :: i

   do i = 1, size(data)
      data(i) = sin(0.19_dp * real(i, dp)) + 0.25_dp * cos(0.73_dp * real(i, dp)) + &
         0.05_dp * real(mod(11 * i, 7) - 3, dp)
   end do
   estimate = bayes_thr(data, filter_number=1.0_dp, family="DaubExPhase", j0=1)
   call require(size(estimate) == size(data), "BAYES.THR reconstruction length")
   call require(all(ieee_is_finite(estimate)), "BAYES.THR finite reconstruction")
   call require(maxval(abs(estimate - data)) > 1.0e-6_dp, "BAYES.THR should shrink the noisy deterministic series")
   call require(maxval(abs(estimate(1:4) - 0.36766441908359937_dp)) < 2.0e-12_dp, &
      "BAYES.THR deterministic Haar regression values")
   call require(maxval(abs(estimate(5:8) - 0.90340515953677847_dp)) < 2.0e-12_dp, &
      "BAYES.THR deterministic second Haar block")

   zeros = 0.0_dp
   zero_estimate = bayes_thr(zeros, filter_number=1.0_dp, family="DaubExPhase", j0=1)
   call require(size(zero_estimate) == size(zeros), "BAYES.THR zero-data length")
   call require(maxval(abs(zero_estimate)) <= tiny(1.0_dp), "BAYES.THR zero-data fixed point")

   print *, "test_bayes: PASS"

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Condition that must hold for the deterministic Bayes threshold test.
      character(len=*), intent(in) :: message !! Failure explanation printed before terminating the test.

      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine require

end program test_bayes
