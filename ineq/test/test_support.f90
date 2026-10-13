! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module test_support
   use ineq_api, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   implicit none
   integer :: checked = 0, failed = 0
contains
   subroutine check(actual, expected, label, tolerance)
      real(dp), intent(in) :: actual !! Computed scalar being tested.
      real(dp), intent(in) :: expected !! Frozen result from R or a mathematical identity.
      character(*), intent(in) :: label !! Description identifying a failed check.
      real(dp), intent(in), optional :: tolerance !! Relative and absolute tolerance, default 2e-12.
      real(dp) :: tol
      logical :: same
      checked = checked + 1
      tol = 2.0e-12_dp
      if (present(tolerance)) tol = tolerance
      if (ieee_is_nan(expected)) then
         same = ieee_is_nan(actual)
      else if (.not. ieee_is_finite(expected)) then
         same = actual == expected
      else
         same = ieee_is_finite(actual)
         if (same) same = abs(actual - expected) <= tol * max(1.0_dp, abs(expected))
      end if
      if (.not. same) then
         print *, 'FAIL: ', label
         print *, 'actual: ', actual, ' expected: ', expected
         failed = failed + 1
      end if
   end subroutine check

   subroutine check_array(actual, expected, label)
      real(dp), intent(in) :: actual(:) !! Computed vector, required to have the expected size.
      real(dp), intent(in) :: expected(:) !! R reference or independently derived expected vector.
      character(*), intent(in) :: label !! Description identifying the vector being checked.
      integer :: i
      call require(size(actual) == size(expected), label // ' length')
      do i = 1, size(expected)
         call check(actual(i), expected(i), label)
      end do
   end subroutine check_array

   subroutine require(condition, label)
      logical, intent(in) :: condition !! Assertion that must evaluate to true.
      character(*), intent(in) :: label !! Explanation reported if the assertion fails.
      checked = checked + 1
      if (.not. condition) then
         print *, 'FAIL: ', label
         error stop 1
      end if
   end subroutine require
end module test_support
