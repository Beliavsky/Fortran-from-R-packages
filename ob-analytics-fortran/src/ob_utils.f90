! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of obAnalytics.
! Copyright (C) 2015,2016 Philip Stubbings.
module ob_utils
   use ob_kinds, only : dp, i8
   implicit none
   private
   public :: vector_diff, reverse_matrix, normalize, interval_sum_breaks, interval_vwap
   public :: weighted_average, interval_price_level_gaps, round_digits, sort_indices_i8
   public :: sort_indices_real_i8, sort_indices_price_time
   public :: unique_i8, contains_i8, find_i8, append_i8

contains

   pure function vector_diff(v) result(d)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: d(:)
      integer :: n
      n = size(v)
      allocate(d(n))
      if (n == 0) return
      d(1) = 0.0_dp
      if (n > 1) d(2:n) = v(2:n) - v(1:n-1)
   end function vector_diff


   pure function reverse_matrix(m) result(reversed)
      real(dp), intent(in) :: m(:,:)
      real(dp), allocatable :: reversed(:,:)
      integer :: i, n
      n = size(m,1)
      allocate(reversed(n,size(m,2)))
      do i = 1, n
         reversed(i,:) = m(n-i+1,:)
      end do
   end function reverse_matrix

   pure function normalize(v, min_value, max_value) result(z)
      real(dp), intent(in) :: v(:)
      real(dp), intent(in), optional :: min_value, max_value
      real(dp), allocatable :: z(:)
      real(dp) :: lo, hi
      allocate(z(size(v)))
      if (size(v) == 0) return
      lo = minval(v)
      hi = maxval(v)
      if (present(min_value)) lo = min_value
      if (present(max_value)) hi = max_value
      if (hi == lo) then
         z = 0.0_dp
      else
         z = (v - lo)/(hi - lo)
      end if
   end function normalize

   pure function interval_sum_breaks(v, breaks) result(sums)
      real(dp), intent(in) :: v(:)
      integer, intent(in) :: breaks(:)
      real(dp), allocatable :: sums(:)
      integer :: i, first, last
      allocate(sums(size(breaks)))
      sums = 0.0_dp
      first = 1
      do i = 1, size(breaks)
         last = min(max(breaks(i), 0), size(v))
         if (last >= first) sums(i) = sum(v(first:last))
         first = last + 1
      end do
   end function interval_sum_breaks


   pure function interval_price_level_gaps(volume, breaks) result(gaps)
      real(dp), intent(in) :: volume(:)
      integer, intent(in) :: breaks(:)
      real(dp), allocatable :: gaps(:)
      gaps = interval_sum_breaks(merge(1.0_dp,0.0_dp,volume==0.0_dp),breaks)
   end function interval_price_level_gaps

   pure function weighted_average(value, weight) result(answer)
      real(dp), intent(in) :: value(:), weight(:)
      real(dp) :: answer
      if (size(value) /= size(weight) .or. size(value) == 0 .or. sum(weight) == 0.0_dp) then
         answer = 0.0_dp
      else
         answer = dot_product(value, weight)/sum(weight)
      end if
   end function weighted_average

   pure function interval_vwap(price, volume, breaks) result(answer)
      real(dp), intent(in) :: price(:), volume(:)
      integer, intent(in) :: breaks(:)
      real(dp), allocatable :: answer(:)
      real(dp), allocatable :: numer(:), denom(:)
      integer :: i
      numer = interval_sum_breaks(price*volume, breaks)
      denom = interval_sum_breaks(volume, breaks)
      allocate(answer(size(breaks)))
      answer = 0.0_dp
      do i = 1, size(answer)
         if (denom(i) /= 0.0_dp) answer(i) = numer(i)/denom(i)
      end do
   end function interval_vwap

   pure elemental function round_digits(x, digits) result(y)
      real(dp), intent(in) :: x
      integer, intent(in) :: digits
      real(dp) :: y, scale
      scale = 10.0_dp**digits
      y = anint(x*scale)/scale
   end function round_digits

   pure function sort_indices_i8(primary, secondary) result(idx)
      integer(i8), intent(in) :: primary(:)
      integer, intent(in), optional :: secondary(:)
      integer, allocatable :: idx(:)
      integer :: i, j, key
      allocate(idx(size(primary)))
      idx = [(i, i=1, size(primary))]
      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (primary(idx(j)) < primary(key)) exit
            if (primary(idx(j)) == primary(key)) then
               if (.not. present(secondary)) exit
               if (secondary(idx(j)) <= secondary(key)) exit
            end if
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end function sort_indices_i8

   pure function sort_indices_real_i8(primary, secondary) result(idx)
      real(dp), intent(in) :: primary(:)
      integer(i8), intent(in) :: secondary(:)
      integer, allocatable :: idx(:)
      integer :: i, j, key
      allocate(idx(size(primary)))
      idx = [(i, i=1, size(primary))]
      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (primary(idx(j)) < primary(key)) exit
            if (primary(idx(j)) == primary(key) .and. secondary(idx(j)) <= secondary(key)) exit
            idx(j+1) = idx(j)
            j = j - 1
         end do
         idx(j+1) = key
      end do
   end function sort_indices_real_i8

   pure function sort_indices_price_time(price, timestamp) result(idx)
      real(dp), intent(in) :: price(:)
      integer(i8), intent(in) :: timestamp(:)
      integer, allocatable :: idx(:)
      idx = sort_indices_real_i8(price, timestamp)
   end function sort_indices_price_time

   pure function unique_i8(v) result(u)
      integer(i8), intent(in) :: v(:)
      integer(i8), allocatable :: u(:)
      integer :: i, n
      allocate(u(size(v)))
      n = 0
      do i = 1, size(v)
         if (.not. any(u(1:n) == v(i))) then
            n = n + 1
            u(n) = v(i)
         end if
      end do
      if (n < size(u)) u = u(1:n)
   end function unique_i8

   pure logical function contains_i8(v, value)
      integer(i8), intent(in) :: v(:), value
      contains_i8 = any(v == value)
   end function contains_i8

   pure integer function find_i8(v, value)
      integer(i8), intent(in) :: v(:), value
      integer :: i
      find_i8 = 0
      do i = 1, size(v)
         if (v(i) == value) then
            find_i8 = i
            return
         end if
      end do
   end function find_i8

   subroutine append_i8(v, value)
      integer(i8), allocatable, intent(inout) :: v(:)
      integer(i8), intent(in) :: value
      integer(i8), allocatable :: tmp(:)
      integer :: n
      if (.not. allocated(v)) then
         allocate(v(1)); v(1) = value; return
      end if
      n = size(v)
      allocate(tmp(n+1))
      if (n > 0) tmp(1:n) = v
      tmp(n+1) = value
      call move_alloc(tmp, v)
   end subroutine append_i8

end module ob_utils
