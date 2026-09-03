! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_utils
   use r_kinds, only : dp
   implicit none
   private
   public :: div_zero
   public :: pair_index_upper
   public :: qform
   public :: trace_matrix

contains

   pure elemental real(dp) function div_zero(x, y, tolerance) result(value)
      real(dp), intent(in) :: x !! Numerator of the requested ratio.
      real(dp), intent(in) :: y !! Denominator of the requested ratio.
      real(dp), intent(in), optional :: tolerance !! Threshold below which two near-zero operands map to one.
      real(dp) :: tol

      tol = 1.0e-14_dp
      if (present(tolerance)) tol = tolerance
      if (abs(x) < tol .and. abs(y) < tol) then
         value = 1.0_dp
      else
         value = x / y
      end if
   end function div_zero

   pure elemental integer function pair_index_upper(i, j, n) result(index_value)
      integer, intent(in) :: i !! One-based row index in a symmetric `n` by `n` matrix.
      integer, intent(in) :: j !! One-based column index in a symmetric `n` by `n` matrix.
      integer, intent(in) :: n !! Order of the symmetric matrix; must be positive.
      integer :: lower_index
      integer :: upper_index

      lower_index = min(i, j)
      upper_index = max(i, j)
      index_value = ((lower_index - 1) * (2 * n - lower_index)) / 2 + upper_index
   end function pair_index_upper

   pure real(dp) function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal is summed.
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         value = value + a(i, i)
      end do
   end function trace_matrix

   pure real(dp) function qform(x, a) result(value)
      real(dp), intent(in) :: x(:) !! Vector defining the quadratic form.
      real(dp), intent(in) :: a(:, :) !! Square matrix conformable with `x`.

      value = dot_product(x, matmul(a, x))
   end function qform

end module pbkrtest_utils
