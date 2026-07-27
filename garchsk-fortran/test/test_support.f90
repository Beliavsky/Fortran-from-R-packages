! SPDX-License-Identifier: GPL-2.0-or-later
module test_support
   use garchsk_kinds, only : dp
   implicit none
   private
   public :: assert_close, assert_true
contains
   subroutine assert_close(actual, expected, atol, rtol)
      real(dp), intent(in) :: actual, expected
      real(dp), intent(in), optional :: atol, rtol
      real(dp) :: abs_tol, rel_tol, error
      abs_tol = 1.0e-12_dp
      rel_tol = 1.0e-10_dp
      if (present(atol)) abs_tol = atol
      if (present(rtol)) rel_tol = rtol
      error = abs(actual - expected)
      if (error > abs_tol + rel_tol * abs(expected)) then
         write(*, '(a,3(1x,es24.16))') 'mismatch:', actual, expected, error
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') trim(message)
         error stop 1
      end if
   end subroutine assert_true
end module test_support
