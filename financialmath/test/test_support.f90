! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module test_support
   use financialmath_kinds, only : dp
   implicit none
   private
   public :: assert_close, assert_true
contains
   subroutine assert_close(actual, expected, atol, rtol)
      real(dp), intent(in) :: actual, expected
      real(dp), intent(in), optional :: atol, rtol
      real(dp) :: a, r, difference
      a = 1.0e-11_dp
      r = 1.0e-10_dp
      if (present(atol)) a = atol
      if (present(rtol)) r = rtol
      difference = abs(actual-expected)
      if (difference > a+r*max(abs(actual), abs(expected))) then
         write(*,'(a,3(1x,es24.16))') 'mismatch:', actual, expected, difference
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') trim(message)
         error stop 1
      end if
   end subroutine assert_true
end module test_support
