! SPDX-License-Identifier: MIT
module r_sorting
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_ordering, only : r_order, r_sort_values_in_place
   use r_quantiles, only : core_quantile_type7 => r_quantile_type7
   use r_quantiles, only : core_weighted_quantile_ecdf => r_weighted_quantile_ecdf
   use r_status, only : r_invalid_input, r_ok
   implicit none
   private

   public :: r_average_ranks, r_quantile_type7, r_sort, r_weighted_quantile_ecdf

contains

   pure subroutine r_sort(x, values, na_rm, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: values(:)
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
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

      call r_sort_values_in_place(values)
   end subroutine r_sort

   pure real(dp) function r_quantile_type7(x, probability, na_rm) result(value)
      real(dp), intent(in) :: x(:), probability
      logical, intent(in), optional :: na_rm
      value = core_quantile_type7(x, probability, na_rm)
   end function r_quantile_type7

   pure real(dp) function r_weighted_quantile_ecdf(x, weights, probability, na_rm, finite_only) result(value)
      real(dp), intent(in) :: x(:), weights(:), probability
      logical, intent(in), optional :: na_rm, finite_only
      value = core_weighted_quantile_ecdf(x, weights, probability, na_rm, finite_only)
   end function r_weighted_quantile_ecdf

   pure subroutine r_average_ranks(x, ranks, tie_tolerance, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: ranks(:)
      real(dp), intent(in), optional :: tie_tolerance
      integer, intent(out), optional :: status
      integer, allocatable :: order(:)
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
      allocate(ranks(n))
      if (n == 0) return
      call r_order(x, order)
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

end module r_sorting
