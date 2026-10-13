! SPDX-License-Identifier: GPL-2.0-or-later
module irr_common
   use r_kinds, only: dp
   use r_ordering, only: r_order, r_sort_values_in_place
   use r_descriptive, only: r_mean, r_variance, r_correlation
   use r_distributions, only: r_pnorm, r_pchisq, r_pf, r_qf
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   private
   integer, parameter, public :: irr_ok = 0, irr_invalid = 1, irr_degenerate = 2, irr_resource = 3
   integer, parameter, public :: max_categories = 256, max_raters = 1000, max_cells = 10000000
   type, public :: irr_result
      real(dp) :: value = 0.0_dp, statistic = 0.0_dp, p_value = 0.0_dp
      real(dp) :: df1 = 0.0_dp, df2 = 0.0_dp, lbound = 0.0_dp, ubound = 0.0_dp
      integer :: status = irr_invalid, subjects = 0, raters = 0, dropped = 0
      real(dp), allocatable :: detail(:, :), levels(:), coincidence(:, :)
   end type
   public :: dp, init_result, complete_ratings, encode_categories, contingency
   public :: average_ranks, mean_squares, finish_z, nan_value, divide, valid_table
contains
   pure function nan_value() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function

   pure subroutine init_result(res)
      type(irr_result), intent(out) :: res !! Result initialized to invalid status and quiet NaN numeric outputs.
      res%value = nan_value()
      res%statistic = nan_value()
      res%p_value = nan_value()
      res%df1 = nan_value()
      res%df2 = nan_value()
      res%lbound = nan_value()
      res%ubound = nan_value()
   end subroutine

   pure subroutine complete_ratings(ratings, x, res, min_subjects)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater matrix; any NaN removes the whole subject.
      real(dp), allocatable, intent(out) :: x(:, :) !! Finite complete-case matrix, or an empty matrix on failure.
      type(irr_result), intent(out) :: res !! Initialized status and counts of retained subjects and raters.
      integer, intent(in), optional :: min_subjects !! Minimum retained subjects; defaults to one.
      integer :: i, j, n, need
      call init_result(res)
      allocate(x(0, 0))
      res%raters = size(ratings, 2)
      if (res%raters < 2) return
      if (res%raters > max_raters .or. size(ratings, kind = 8) > int(max_cells, 8)) then
         res%status = irr_resource
         return
      end if
      n = 0
      do i = 1, size(ratings, 1)
         if (any(ieee_is_nan(ratings(i, :)))) cycle
         if (any(.not. ieee_is_finite(ratings(i, :)))) return
         n = n + 1
      end do
      res%subjects = n
      need = 1
      if (present(min_subjects)) need = min_subjects
      if (n < need) return
      deallocate(x)
      allocate(x(n, res%raters))
      j = 0
      do i = 1, size(ratings, 1)
         if (any(ieee_is_nan(ratings(i, :)))) cycle
         j = j + 1
         x(j, :) = ratings(i, :)
      end do
      res%status = irr_ok
   end subroutine

   pure subroutine encode_categories(x, codes, levels, status)
      real(dp), intent(in) :: x(:, :) !! Finite numeric category labels, sorted in increasing numeric order.
      integer, allocatable, intent(out) :: codes(:, :) !! One-based category codes with the shape of x.
      real(dp), allocatable, intent(out) :: levels(:) !! Distinct numeric category labels in ascending order.
      integer, intent(out) :: status !! irr_ok or irr_resource when more than 256 categories are present.
      real(dp), allocatable :: work(:)
      integer :: i, j, n, k, lo, hi, mid
      status = irr_resource
      allocate(codes(0, 0), levels(0))
      if (size(x, kind = 8) > int(max_cells, 8)) return
      if (any(.not. ieee_is_finite(x))) then
         status = irr_invalid
         return
      end if
      work = reshape(x, [size(x)])
      call r_sort_values_in_place(work)
      n = 0
      do i = 1, size(work)
         if (i > 1) then
            if (work(i) == work(i - 1)) cycle
         end if
         n = n + 1
         if (n > max_categories) return
      end do
      deallocate(levels, codes)
      allocate(levels(n), codes(size(x, 1), size(x, 2)))
      k = 0
      do i = 1, size(work)
         if (i > 1) then
            if (work(i) == work(i - 1)) cycle
         end if
         k = k + 1
         levels(k) = work(i)
      end do
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            lo = 1
            hi = n
            do while (lo < hi)
               mid = lo + (hi - lo) / 2
               if (levels(mid) < x(i, j)) then
                  lo = mid + 1
               else
                  hi = mid
               end if
            end do
            codes(i, j) = lo
         end do
      end do
      status = irr_ok
   end subroutine

   pure subroutine contingency(ratings, tab, res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-two-rater category matrix with NaN listwise deletion.
      real(dp), allocatable, intent(out) :: tab(:, :) !! Square count table using a shared sorted category universe.
      type(irr_result), intent(out) :: res !! Status, counts, and category labels.
      real(dp), allocatable :: x(:, :)
      integer, allocatable :: codes(:, :)
      integer :: i, nc
      allocate(tab(0, 0))
      call complete_ratings(ratings, x, res)
      if (res%status /= irr_ok) return
      if (res%raters /= 2) then
         res%status = irr_invalid
         return
      end if
      call encode_categories(x, codes, res%levels, res%status)
      if (res%status /= irr_ok) return
      nc = size(res%levels)
      deallocate(tab)
      allocate(tab(nc, nc))
      tab = 0.0_dp
      do i = 1, size(x, 1)
         tab(codes(i, 1), codes(i, 2)) = tab(codes(i, 1), codes(i, 2)) + 1.0_dp
      end do
   end subroutine

   pure subroutine average_ranks(x, ranks, ties)
      real(dp), intent(in) :: x(:) !! Finite observations to rank in ascending order.
      real(dp), intent(out) :: ranks(:) !! Average ranks for ties; must have the same size as x.
      real(dp), intent(out) :: ties !! Sum of t**3 - t over tie groups of size t.
      integer, allocatable :: order(:)
      integer :: i, j, t
      call r_order(x, order)
      ties = 0.0_dp
      i = 1
      do while (i <= size(x))
         j = i
         do while (j < size(x))
            if (x(order(j + 1)) /= x(order(i))) exit
            j = j + 1
         end do
         ranks(order(i:j)) = real(i + j, dp) / 2.0_dp
         t = j - i + 1
         ties = ties + real(t, dp) ** 3 - real(t, dp)
         i = j + 1
      end do
   end subroutine

   pure subroutine mean_squares(x, msr, msw, msc, mse)
      real(dp), intent(in) :: x(:, :) !! Complete finite subject-by-rater matrix, both dimensions at least two.
      real(dp), intent(out) :: msr !! Subject mean square.
      real(dp), intent(out) :: msw !! Pooled within-subject mean square.
      real(dp), intent(out) :: msc !! Rater mean square.
      real(dp), intent(out) :: mse !! Two-way residual mean square.
      real(dp) :: rows(size(x, 1)), cols(size(x, 2)), grand, total
      integer :: n, r, i
      n = size(x, 1)
      r = size(x, 2)
      rows = sum(x, dim = 2) / real(r, dp)
      cols = sum(x, dim = 1) / real(n, dp)
      grand = r_mean(rows)
      total = sum((x - grand) ** 2)
      msr = sum((rows - grand) ** 2) * real(r, dp) / real(n - 1, dp)
      msc = sum((cols - grand) ** 2) * real(n, dp) / real(r - 1, dp)
      msw = 0.0_dp
      mse = 0.0_dp
      do i = 1, n
         msw = msw + sum((x(i, :) - rows(i)) ** 2)
         mse = mse + sum((x(i, :) - rows(i) - cols + grand) ** 2)
      end do
      msw = msw / real(n * (r - 1), dp)
      mse = mse / real((n - 1) * (r - 1), dp)
   end subroutine

   pure elemental function divide(a, b) result(x)
      real(dp), intent(in) :: a !! Finite numerator of a statistic or coefficient.
      real(dp), intent(in) :: b !! Denominator; zero yields signed infinity, or NaN when a is also zero.
      real(dp) :: x
      if (b == 0.0_dp) then
         if (a == 0.0_dp) then
            x = nan_value()
         else
            x = sign(ieee_value(0.0_dp, ieee_positive_inf), a)
         end if
      else
         x = a / b
      end if
   end function

   pure subroutine finish_z(res, variance)
      type(irr_result), intent(inout) :: res !! Coefficient result to augment with its normal statistic and two-sided p-value.
      real(dp), intent(in) :: variance !! Variance of the coefficient; negative or NaN marks degenerate inference.
      if (.not. ieee_is_finite(variance) .or. variance < 0.0_dp .or. ieee_is_nan(res%value)) then
         res%status = irr_degenerate
         return
      end if
      res%statistic = divide(res%value, sqrt(variance))
      res%p_value = 2.0_dp * r_pnorm(abs(res%statistic), lower_tail = .false.)
      if (ieee_is_nan(res%statistic)) res%status = irr_degenerate
   end subroutine

   pure function valid_table(tab) result(ok)
      real(dp), intent(in) :: tab(:, :) !! Nonnegative finite square contingency counts, allowed to be fractional.
      logical :: ok
      ok = .false.
      if (size(tab, 1) < 1 .or. size(tab, 1) /= size(tab, 2)) return
      if (size(tab, 1) > max_categories) return
      if (any(.not. ieee_is_finite(tab))) return
      if (any(tab < 0.0_dp)) return
      if (sum(tab) <= 0.0_dp .or. sum(tab) > real(huge(0), dp)) return
      ok = .true.
   end function
end module
