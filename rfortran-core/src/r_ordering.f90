! SPDX-License-Identifier: MIT
module r_ordering
   use r_kinds, only : dp
   implicit none
   private

   public :: r_order, r_sort_values_in_place

contains

   pure subroutine r_order(x, order)
      real(dp), intent(in) :: x(:)
      integer, allocatable, intent(out) :: order(:)
      integer, allocatable :: work(:)
      integer :: i

      allocate(order(size(x)))
      if (size(x) == 0) return
      order = [(i, i=1,size(x))]
      if (size(x) == 1) return
      allocate(work(size(x)))
      call merge_sort_indices(x, order, work, 1, size(x))
   end subroutine r_order

   pure subroutine r_sort_values_in_place(values)
      real(dp), intent(inout) :: values(:)
      real(dp), allocatable :: work(:)

      if (size(values) < 2) return
      allocate(work(size(values)))
      call merge_sort_values(values, work, 1, size(values))
   end subroutine r_sort_values_in_place

   pure recursive subroutine merge_sort_values(values, work, left, right)
      real(dp), intent(inout) :: values(:), work(:)
      integer, intent(in) :: left, right
      integer :: i, j, k, middle

      if (left >= right) return
      middle = (left + right)/2
      call merge_sort_values(values, work, left, middle)
      call merge_sort_values(values, work, middle + 1, right)
      i = left
      j = middle + 1
      do k = left, right
         if (i > middle) then
            work(k) = values(j)
            j = j + 1
         else if (j > right) then
            work(k) = values(i)
            i = i + 1
         else if (values(i) <= values(j)) then
            work(k) = values(i)
            i = i + 1
         else
            work(k) = values(j)
            j = j + 1
         end if
      end do
      values(left:right) = work(left:right)
   end subroutine merge_sort_values

   pure recursive subroutine merge_sort_indices(x, order, work, left, right)
      real(dp), intent(in) :: x(:)
      integer, intent(inout) :: order(:), work(:)
      integer, intent(in) :: left, right
      integer :: i, j, k, middle

      if (left >= right) return
      middle = (left + right)/2
      call merge_sort_indices(x, order, work, left, middle)
      call merge_sort_indices(x, order, work, middle + 1, right)
      i = left
      j = middle + 1
      do k = left, right
         if (i > middle) then
            work(k) = order(j)
            j = j + 1
         else if (j > right) then
            work(k) = order(i)
            i = i + 1
         else if (x(order(i)) <= x(order(j))) then
            work(k) = order(i)
            i = i + 1
         else
            work(k) = order(j)
            j = j + 1
         end if
      end do
      order(left:right) = work(left:right)
   end subroutine merge_sort_indices

end module r_ordering
