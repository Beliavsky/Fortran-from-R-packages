! SPDX-License-Identifier: GPL-2.0-or-later
! randtests' permut is attributed upstream to gtools::permutations.
module randtests_permutations
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only: dp
   implicit none
   private
   integer, parameter, public :: max_permutation_entries = 1000000
   public :: permut, permut_apply, randtests_aux, permutation_callback
   abstract interface
      pure function permutation_callback(x) result(value)
         import dp
         real(dp), intent(in) :: x(:) !! One ordered selection of m elements from the source vector.
         real(dp) :: value
      end function permutation_callback
   end interface
contains
   pure subroutine permut(x, permutations, status, m)
      real(dp), intent(in) :: x(:) !! Finite source elements; repeated values represent distinct source positions.
      real(dp), allocatable, intent(out) :: permutations(:, :) !! Matrix with one ordered selection in each row.
      integer, intent(out) :: status !! Zero success, one invalid input, two resource limit exceeded.
      integer, intent(in), optional :: m !! Selection length from 0 through size(x); default size(x).
      integer :: length, rows, i, row, cells
      logical, allocatable :: used(:)
      real(dp), allocatable :: current(:)
      length = size(x)
      if (present(m)) length = m
      status = 1
      if (length < 0 .or. length > size(x) .or. any(.not. ieee_is_finite(x))) then
         allocate (permutations(0, 0))
         return
      end if
      rows = 1
      cells = max(1, length)
      do i = 0, length - 1
         if (rows > max_permutation_entries / cells / (size(x) - i)) then
            allocate (permutations(0, 0))
            status = 2
            return
         end if
         rows = rows * (size(x) - i)
      end do
      allocate (permutations(rows, length), current(length), used(size(x)))
      used = .false.
      row = 0
      call enumerate_permutations(x, 1, current, used, permutations, row)
      status = 0
   end subroutine permut

   pure subroutine permut_apply(x, fun, values, status, m)
      real(dp), intent(in) :: x(:) !! Finite source elements for ordered selections without replacement.
      procedure(permutation_callback), pointer, intent(in) :: fun !! Pure scalar-real callback, applied once to every permutation.
      real(dp), allocatable, intent(out) :: values(:) !! Callback values in the same row order as permut.
      integer, intent(out) :: status !! Zero success, one invalid input, two resource limit exceeded.
      integer, intent(in), optional :: m !! Selection length in [0,size(x)]; default size(x).
      real(dp), allocatable :: permutations(:, :)
      integer :: i
      if (.not. associated(fun)) then
         allocate (values(0))
         status = 1
         return
      end if
      call permut(x, permutations, status, m)
      if (status /= 0) then
         allocate (values(0))
         return
      end if
      allocate (values(size(permutations, 1)))
      do i = 1, size(values)
         values(i) = fun(permutations(i, :))
      end do
   end subroutine permut_apply

   pure recursive subroutine enumerate_permutations(x, depth, current, used, permutations, row)
      real(dp), intent(in) :: x(:) !! Original finite source vector, with positions treated as distinct.
      integer, intent(in) :: depth !! Position being filled in the current ordered selection.
      real(dp), intent(inout) :: current(:) !! Working selection of the requested length.
      logical, intent(inout) :: used(:) !! Source-position mask for the current recursion path.
      real(dp), intent(inout) :: permutations(:, :) !! Preallocated output matrix, one selection per row.
      integer, intent(inout) :: row !! Number of output rows filled so far.
      integer :: i
      if (depth > size(current)) then
         row = row + 1
         permutations(row, :) = current
         return
      end if
      do i = 1, size(x)
         if (used(i)) cycle
         used(i) = .true.
         current(depth) = x(i)
         call enumerate_permutations(x, depth + 1, current, used, permutations, row)
         used(i) = .false.
      end do
   end subroutine enumerate_permutations

   pure real(dp) function randtests_aux(x, method) result(value)
      real(dp), intent(in) :: x(:) !! Finite sequence; empty and singleton sequences have statistic zero.
      character(len=*), intent(in) :: method !! Select "bartels" squared differences or "difference.sign" increases.
      integer :: i
      real(dp) :: scale, difference, scaled_sum, root_max
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (any(.not. ieee_is_finite(x))) return
      select case (method)
      case ('difference.sign')
         value = real(count(x(2:) > x(:size(x) - 1)), dp)
      case ('bartels')
         value = 0.0_dp
         if (size(x) < 2) return
         ! Subtract before scaling to preserve differences in large-offset data.
         ! Same-sign subtraction cannot overflow; an opposite-sign difference
         ! exceeding sqrt(huge) already implies an unrepresentable square.
         root_max = sqrt(huge(1.0_dp))
         scale = 0.0_dp
         scaled_sum = 1.0_dp
         do i = 2, size(x)
            if ((x(i) > 0.0_dp) .neqv. (x(i - 1) > 0.0_dp)) then
               if (max(abs(x(i)), abs(x(i - 1))) > root_max) then
                  value = ieee_value(0.0_dp, ieee_quiet_nan)
                  return
               end if
            end if
            difference = abs(x(i) - x(i - 1))
            if (difference > root_max) then
               value = ieee_value(0.0_dp, ieee_quiet_nan)
               return
            end if
            if (difference == 0.0_dp) cycle
            if (scale < difference) then
               scaled_sum = 1.0_dp + scaled_sum * (scale / difference) ** 2
               scale = difference
            else
               scaled_sum = scaled_sum + (difference / scale) ** 2
            end if
         end do
         if (scale == 0.0_dp) return
         if (scale > root_max / sqrt(scaled_sum)) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            value = (scale * sqrt(scaled_sum)) ** 2
         end if
      end select
   end function randtests_aux
end module randtests_permutations
