! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_utils
   use r_kinds, only : dp
   implicit none
   private

   public :: stable_order_real
   public :: unique_sorted_real
   public :: identity_matrix

contains

   pure subroutine stable_order_real(x, order)
      real(dp), intent(in) :: x(:) !! Values whose stable ascending order is requested.
      integer, intent(out) :: order(:) !! Permutation indices, with size equal to `size(x)`.

      integer :: i
      integer :: j
      integer :: key

      do i = 1, size(x)
         order(i) = i
      end do
      do i = 2, size(x)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (x(order(j)) <= x(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine stable_order_real

   pure subroutine unique_sorted_real(x, unique_x)
      real(dp), intent(in) :: x(:) !! Values assumed sorted in nondecreasing order.
      real(dp), allocatable, intent(out) :: unique_x(:) !! Distinct values in ascending order.

      real(dp), allocatable :: work(:)
      integer :: count
      integer :: i

      if (size(x) == 0) then
         allocate(unique_x(0))
         return
      end if

      allocate(work(size(x)))
      count = 1
      work(1) = x(1)
      do i = 2, size(x)
         if (x(i) /= work(count)) then
            count = count + 1
            work(count) = x(i)
         end if
      end do
      allocate(unique_x(count))
      unique_x = work(1:count)
   end subroutine unique_sorted_real

   pure subroutine identity_matrix(a)
      real(dp), intent(out) :: a(:, :) !! Square matrix replaced by the identity matrix.

      integer :: i

      a = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         a(i, i) = 1.0_dp
      end do
   end subroutine identity_matrix

end module cmprsk_utils
