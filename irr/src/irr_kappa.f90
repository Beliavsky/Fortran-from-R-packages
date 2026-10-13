! SPDX-License-Identifier: GPL-2.0-or-later
module irr_kappa
   use irr_common
   use r_ordering, only: r_sort_values_in_place
   use r_distributions, only: r_pnorm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: kappa2, kappa2_table, kappam_fleiss, kappam_light
contains
   pure function kappa2(ratings, weight, costs) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-two-rater numeric categories; incomplete rows are omitted.
      character(*), intent(in), optional :: weight !! 'unweighted' (default), 'equal', or 'squared'; ignored when costs are
      !! supplied.
      real(dp), intent(in), optional :: costs(:) !! One cost per category separation; normalized exactly as the R numeric weight
      !! vector.
      type(irr_result) :: res, base
      real(dp), allocatable :: tab(:, :)
      call contingency(ratings, tab, base)
      res = base
      if (base%status /= irr_ok) return
      res = kappa2_table(tab, weight, costs)
      res%levels = base%levels
   end function

   pure function kappa2_table(tab, weight, costs) result(res)
      real(dp), intent(in) :: tab(:, :) !! Nonnegative square contingency counts with identical ordered row/column category labels.
      character(*), intent(in), optional :: weight !! 'unweighted' (default), 'equal', or 'squared' agreement weights.
      real(dp), intent(in), optional :: costs(:) !! Numeric R-style costs, length equal to table order, with nonzero range.
      type(irr_result) :: res
      real(dp), allocatable :: w(:), wvec(:), wt(:, :), p(:), q(:), e(:, :), wi(:), wj(:)
      real(dp) :: n, observed, chance, variance
      integer :: nc, i, j
      character(:), allocatable :: kind
      call init_result(res)
      if (.not. valid_table(tab)) return
      nc = size(tab, 1)
      n = sum(tab)
      res%subjects = nint(n)
      res%raters = 2
      res%status = irr_degenerate
      if (nc < 2) return
      kind = 'unweighted'
      if (present(weight)) kind = weight
      allocate(w(nc), wvec(2 * nc - 1), wt(nc, nc), e(nc, nc), wi(nc), wj(nc))
      if (present(costs)) then
         res%status = irr_invalid
         if (size(costs) /= nc) return
         if (any(.not. ieee_is_finite(costs))) return
         if (maxval(costs) == minval(costs)) return
         w = 1.0_dp - (costs - minval(costs)) / (maxval(costs) - minval(costs))
      else
         select case (kind)
         case ('unweighted')
            w = 0.0_dp
            w(1) = 1.0_dp
         case ('equal')
            w = [(1.0_dp - real(i - 1, dp) / real(nc - 1, dp), i = 1,nc)]
         case ('squared')
            w = [(1.0_dp - (real(i - 1, dp) / real(nc - 1, dp)) ** 2, i = 1,nc)]
         case default
            res%status = irr_invalid
            return
         end select
      end if
      wvec(1:nc) = w
      call r_sort_values_in_place(wvec(1:nc))
      wvec(nc + 1:) = w(2:)
      do i = 1, nc
         wt(i, :) = wvec(nc - i + 1:2 * nc - i)
      end do
      p = sum(tab, dim = 2) / n
      q = sum(tab, dim = 1) / n
      do j = 1, nc
         e(:, j) = p * q(j)
      end do
      observed = sum(tab * wt) / n
      chance = sum(e * wt)
      res%status = irr_degenerate
      if (chance >= 1.0_dp) return
      res%value = (observed - chance) / (1.0_dp - chance)
      ! These are the exact R apply/recycling orientations, also for asymmetric custom weights.
      do i = 1, nc
         wi(i) = sum(q * wt(:, i))
         wj(i) = sum(p * wt(i, :))
      end do
      variance = 0.0_dp
      do j = 1, nc
         variance = variance + sum(e(:, j) * (wt(:, j) - wi - wj(j)) ** 2)
      end do
      variance = (variance - chance ** 2) / (n * (1.0_dp - chance) ** 2)
      res%status = irr_ok
      call finish_z(res, variance)
   end function

   pure function kappam_fleiss(ratings, exact, detail) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater categorical ratings; every incomplete row is omitted.
      logical, intent(in), optional :: exact !! Use the Conger exact chance agreement correction; default false.
      logical, intent(in), optional :: detail !! Include category-specific unrounded kappa, z, and p-value; default false.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :), counts(:, :), margins(:, :), p(:), q(:), kappa(:)
      integer, allocatable :: codes(:, :)
      integer :: i, j, nc
      real(dp) :: n, r, observed, chance, variance, pq
      logical :: exact_used, detail_used
      call complete_ratings(ratings, x, res)
      if (res%status /= irr_ok) return
      call encode_categories(x, codes, res%levels, res%status)
      if (res%status /= irr_ok) return
      n = real(res%subjects, dp)
      r = real(res%raters, dp)
      nc = size(res%levels)
      allocate(counts(res%subjects, nc), margins(res%raters, nc))
      counts = 0.0_dp
      margins = 0.0_dp
      do j = 1, res%raters
         do i = 1, res%subjects
            counts(i, codes(i, j)) = counts(i, codes(i, j)) + 1.0_dp
            margins(j, codes(i, j)) = margins(j, codes(i, j)) + 1.0_dp / n
         end do
      end do
      observed = (sum(counts ** 2) - n * r) / (n * r * (r - 1.0_dp))
      p = sum(counts, dim = 1) / (n * r)
      q = 1.0_dp - p
      chance = sum(p ** 2)
      exact_used = .false.
      detail_used = .false.
      if (present(exact)) exact_used = exact
      if (present(detail)) detail_used = detail
      if (exact_used) then
         do j = 1, nc
            chance = chance - sum((margins(:, j) - p(j)) ** 2) / (r * (r - 1.0_dp))
         end do
      end if
      if (chance >= 1.0_dp) then
         res%status = irr_degenerate
         return
      end if
      res%value = (observed - chance) / (1.0_dp - chance)
      if (exact_used) return
      pq = sum(p * q)
      variance = 2.0_dp * (pq ** 2 - sum(p * q * (q - p))) / (pq ** 2 * n * r * (r - 1.0_dp))
      call finish_z(res, max(0.0_dp, variance))
      if (detail_used) then
         allocate(res%detail(nc, 3))
         kappa = (sum(counts ** 2, dim = 1) - n * r * p) / (n * r * (r - 1.0_dp) * p)
         res%detail(:, 1) = (kappa - p) / q
         res%detail(:, 2) = res%detail(:, 1) / sqrt(2.0_dp / (n * r * (r - 1.0_dp)))
         res%detail(:, 3) = 2.0_dp * r_pnorm(abs(res%detail(:, 2)), lower_tail = .false.)
      end if
   end function

   pure function kappam_light(ratings) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-rater categorical ratings; complete rows are used for every rater pair.
      type(irr_result) :: res, pair
      real(dp), allocatable :: x(:, :), levels(:), margins(:, :)
      integer, allocatable :: codes(:, :)
      real(dp) :: total, disagreement, log_product, b_scaled, variance
      integer :: i, j, k, nc, pairs
      call complete_ratings(ratings, x, res)
      if (res%status /= irr_ok) return
      call encode_categories(x, codes, levels, res%status)
      if (res%status /= irr_ok) return
      nc = size(levels)
      allocate(margins(res%raters, nc))
      margins = 0.0_dp
      do j = 1, res%raters
         do i = 1, res%subjects
            margins(j, codes(i, j)) = margins(j, codes(i, j)) + 1.0_dp / real(res%subjects, dp)
         end do
      end do
      total = 0.0_dp
      pairs = 0
      log_product = 0.0_dp
      do i = 1, res%raters - 1
         do j = i + 1, res%raters
            pair = kappa2(x(:, [i, j]))
            if (ieee_is_nan(pair%value)) then
               res%status = irr_degenerate
               return
            end if
            total = total + pair%value
            pairs = pairs + 1
            disagreement = 1.0_dp - sum(margins(i, :) * margins(j, :))
            if (disagreement <= 0.0_dp) then
               res%status = irr_degenerate
               return
            end if
            log_product = log_product + log(disagreement)
         end do
      end do
      res%value = total / real(pairs, dp)
      b_scaled = real(pairs, dp) * exp(log_product)
      variance = divide(1.0_dp - b_scaled, real(res%subjects, dp) * b_scaled)
      call finish_z(res, variance)
   end function
end module
