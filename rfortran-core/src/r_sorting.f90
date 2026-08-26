! SPDX-License-Identifier: MIT
module r_sorting
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_status, only : r_invalid_input, r_ok
   implicit none
   private

   public :: r_average_ranks, r_quantile_type7, r_sort

contains

   pure subroutine r_sort(x, values, na_rm, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
      real(dp), allocatable :: work(:)
      integer :: i, n
      logical :: remove_na

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      if (present(status)) status = r_ok
      if (.not. remove_na .and. any(ieee_is_nan(x))) then
         if (present(status)) status = r_invalid_input
         return
      end if

      if (remove_na) then
         n = count(.not. ieee_is_nan(x))
         allocate(values(n))
         n = 0
         do i = 1, size(x)
            if (.not. ieee_is_nan(x(i))) then
               n = n + 1
               values(n) = x(i)
            end if
         end do
      else
         allocate(values(size(x)))
         values = x
      end if
      if (size(values) < 2) return

      allocate(work(size(values)))
      call merge_sort_values(values, work, 1, size(values))
   end subroutine r_sort

   pure real(dp) function r_quantile_type7(x, probability, na_rm) result(value)
      real(dp), intent(in) :: x(:), probability
      logical, intent(in), optional :: na_rm
      real(dp), allocatable :: sorted(:)
      real(dp) :: fraction, h
      integer :: j, local_status, n

      if (probability < 0.0_dp .or. probability > 1.0_dp .or. ieee_is_nan(probability)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      call r_sort(x, sorted, na_rm, local_status)
      if (local_status /= r_ok) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      n = size(sorted)
      if (n == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      h = 1.0_dp + real(n - 1, dp)*probability
      j = int(floor(h))
      fraction = h - real(j, dp)
      if (j >= n) then
         value = sorted(n)
      else
         value = (1.0_dp - fraction)*sorted(j) + fraction*sorted(j + 1)
      end if
   end function r_quantile_type7

   pure subroutine r_average_ranks(x, ranks, tie_tolerance, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: ranks(:)
      real(dp), intent(in), optional :: tie_tolerance
      integer, intent(out), optional :: status
      integer, allocatable :: order(:), work(:)
      integer :: i, left, n, right
      real(dp) :: average_rank, scale, tolerance

      if (present(status)) status = r_ok
      if (any(ieee_is_nan(x))) then
         if (present(status)) status = r_invalid_input
         return
      end if
      tolerance = 0.0_dp
      if (present(tie_tolerance)) tolerance = tie_tolerance
      if (tolerance < 0.0_dp) then
         if (present(status)) status = r_invalid_input
         return
      end if

      n = size(x)
      allocate(ranks(n), order(n), work(n))
      if (n == 0) return
      order = [(i, i=1,n)]
      call merge_sort_indices(x, order, work, 1, n)
      scale = max(1.0_dp, maxval(abs(x)))
      left = 1
      do while (left <= n)
         right = left
         do while (right < n)
            if (abs(x(order(right + 1)) - x(order(left))) > tolerance*scale) exit
            right = right + 1
         end do
         average_rank = 0.5_dp*real(left + right, dp)
         do i = left, right
            ranks(order(i)) = average_rank
         end do
         left = right + 1
      end do
   end subroutine r_average_ranks

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

end module r_sorting
