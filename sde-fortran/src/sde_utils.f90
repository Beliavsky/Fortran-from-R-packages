! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use sde_kinds, only : dp
   implicit none
   private

   public :: arithmetic_mean
   public :: sample_standard_deviation
   public :: linspace
   public :: clamp
   public :: all_finite
   public :: interpolate_missing
   public :: sort_indices

   interface clamp
      module procedure clamp_real
   end interface clamp

contains

   pure function arithmetic_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value

      if (size(x) == 0) error stop "arithmetic_mean: empty input"
      value = sum(x)/real(size(x), dp)
   end function arithmetic_mean

   pure function sample_standard_deviation(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      real(dp) :: mean_value

      if (size(x) < 2) then
         value = 0.0_dp
         return
      end if
      mean_value = arithmetic_mean(x)
      value = sqrt(sum((x-mean_value)**2)/real(size(x)-1, dp))
   end function sample_standard_deviation

   pure function linspace(a, b, n) result(values)
      real(dp), intent(in) :: a
      real(dp), intent(in) :: b
      integer, intent(in) :: n
      real(dp), allocatable :: values(:)
      integer :: i

      if (n <= 0) error stop "linspace: n must be positive"
      allocate(values(n))
      if (n == 1) then
         values(1) = a
      else
         do i = 1, n
            values(i) = a+(b-a)*real(i-1, dp)/real(n-1, dp)
         end do
      end if
   end function linspace

   elemental pure function clamp_real(x, lower, upper) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: lower
      real(dp), intent(in) :: upper
      real(dp) :: value
      value = max(lower, min(upper, x))
   end function clamp_real

   pure function all_finite(x) result(value)
      real(dp), intent(in) :: x(:)
      logical :: value
      integer :: i

      value = .true.
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) then
            value = .false.
            return
         end if
      end do
   end function all_finite

   subroutine interpolate_missing(x)
      real(dp), intent(inout) :: x(:)
      integer :: n, first_valid, last_valid, i, j, left, right
      real(dp) :: fraction

      n = size(x)
      if (n == 0) return
      first_valid = 0
      do i = 1, n
         if (.not. ieee_is_nan(x(i))) then
            first_valid = i
            exit
         end if
      end do
      if (first_valid == 0) error stop "interpolate_missing: all values are NaN"
      x(1:first_valid-1) = x(first_valid)

      last_valid = first_valid
      i = first_valid+1
      do while (i <= n)
         if (.not. ieee_is_nan(x(i))) then
            last_valid = i
            i = i+1
         else
            left = i-1
            right = i
            do while (right <= n .and. ieee_is_nan(x(right)))
               right = right+1
            end do
            if (right > n) then
               x(i:n) = x(left)
               exit
            end if
            do j = i, right-1
               fraction = real(j-left, dp)/real(right-left, dp)
               x(j) = (1.0_dp-fraction)*x(left)+fraction*x(right)
            end do
            last_valid = right
            i = right+1
         end if
      end do
   end subroutine interpolate_missing

   subroutine sort_indices(values, indices)
      real(dp), intent(in) :: values(:)
      integer, intent(out) :: indices(:)
      integer :: i, j, key

      if (size(indices) /= size(values)) error stop "sort_indices: size mismatch"
      do i = 1, size(values)
         indices(i) = i
      end do
      do i = 2, size(values)
         key = indices(i)
         j = i-1
         do while (j >= 1)
            if (values(indices(j)) <= values(key)) exit
            indices(j+1) = indices(j)
            j = j-1
         end do
         indices(j+1) = key
      end do
   end subroutine sort_indices

end module sde_utils
