! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_sort
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: sort_real, sort_real_with_index, median, quantile_type7, weighted_high_median
contains
   recursive subroutine quicksort_pair(x, idx, lo, hi)
      real(dp), intent(inout) :: x(:)
      integer, intent(inout) :: idx(:)
      integer, intent(in) :: lo, hi
      integer :: i, j, it
      real(dp) :: pivot, t
      if (lo >= hi) return
      i = lo
      j = hi
      pivot = x((lo + hi) / 2)
      do
         do while (x(i) < pivot)
            i = i + 1
         end do
         do while (x(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            t = x(i); x(i) = x(j); x(j) = t
            it = idx(i); idx(i) = idx(j); idx(j) = it
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (lo < j) call quicksort_pair(x, idx, lo, j)
      if (i < hi) call quicksort_pair(x, idx, i, hi)
   end subroutine quicksort_pair

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer, allocatable :: idx(:)
      integer :: i
      allocate(idx(size(x)))
      idx = [(i, i=1,size(x))]
      if (size(x) > 1) call quicksort_pair(x, idx, 1, size(x))
   end subroutine sort_real

   subroutine sort_real_with_index(x, idx)
      real(dp), intent(inout) :: x(:)
      integer, intent(out) :: idx(:)
      integer :: i
      if (size(idx) /= size(x)) error stop "sort_real_with_index: size mismatch"
      idx = [(i, i=1,size(x))]
      if (size(x) > 1) call quicksort_pair(x, idx, 1, size(x))
   end subroutine sort_real_with_index

   function median(x) result(m)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      real(dp), allocatable :: y(:)
      integer :: n
      n = size(x)
      if (n == 0) then
         m = 0.0_dp
         return
      end if
      y = x
      call sort_real(y)
      if (mod(n,2) == 1) then
         m = y((n+1)/2)
      else
         m = 0.5_dp * (y(n/2) + y(n/2+1))
      end if
   end function median

   function quantile_type7(x, prob) result(q)
      real(dp), intent(in) :: x(:), prob
      real(dp) :: q, h, g
      real(dp), allocatable :: y(:)
      integer :: n, j
      n = size(x)
      if (n == 0) then
         q = 0.0_dp
         return
      end if
      y = x
      call sort_real(y)
      if (prob <= 0.0_dp) then
         q = y(1)
      else if (prob >= 1.0_dp) then
         q = y(n)
      else
         h = 1.0_dp + real(n-1,dp) * prob
         j = floor(h)
         g = h - real(j,dp)
         if (j >= n) then
            q = y(n)
         else
            q = (1.0_dp-g)*y(j) + g*y(j+1)
         end if
      end if
   end function quantile_type7

   function weighted_high_median(x, weights) result(m)
      real(dp), intent(in) :: x(:), weights(:)
      real(dp) :: m, total, acc
      real(dp), allocatable :: y(:), w(:)
      integer, allocatable :: idx(:)
      integer :: i
      if (size(x) /= size(weights) .or. size(x) == 0) error stop "weighted_high_median: invalid input"
      if (any(weights < 0.0_dp)) error stop "weighted_high_median: negative weight"
      y = x
      allocate(idx(size(x)), w(size(x)))
      call sort_real_with_index(y, idx)
      do i = 1, size(x)
         w(i) = weights(idx(i))
      end do
      total = sum(w)
      if (total <= 0.0_dp) then
         m = y(size(y))
         return
      end if
      acc = 0.0_dp
      do i = 1, size(y)
         acc = acc + w(i)
         if (acc > 0.5_dp*total) then
            m = y(i)
            return
         end if
      end do
      m = y(size(y))
   end function weighted_high_median
end module robustbase_sort
