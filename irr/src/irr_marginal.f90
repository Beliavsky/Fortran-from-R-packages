! SPDX-License-Identifier: GPL-2.0-or-later
module irr_marginal
   use irr_common
   use r_linalg, only: solve_system
   use r_distributions, only: r_pchisq
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private
   public :: bhapkar, bhapkar_table, rater_bias, rater_bias_ratings, stuart_maxwell_mh, stuart_maxwell_ratings
contains
   pure function bhapkar(ratings) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-two-rater category matrix with NaN listwise deletion.
      type(irr_result) :: res, base
      real(dp), allocatable :: tab(:, :)
      call contingency(ratings, tab, base)
      res = base
      if (base%status /= irr_ok) return
      res = bhapkar_table(tab)
      res%levels = base%levels
   end function

   pure function bhapkar_table(tab) result(res)
      real(dp), intent(in) :: tab(:, :) !! Square nonnegative contingency counts in a shared ordered category universe.
      type(irr_result) :: res
      res = homogeneity(tab, .true., .false.)
   end function

   pure function stuart_maxwell_mh(tab, legacy_equal_marginals) result(res)
      real(dp), intent(in) :: tab(:, :) !! Square contingency counts; unlike R, rectangular ratings are never guessed from shape.
      logical, intent(in), optional :: legacy_equal_marginals !! Default true removes equal-marginal categories as R 0.85 does;
      !! false retains discordant categories.
      type(irr_result) :: res
      logical :: legacy
      legacy = .true.
      if (present(legacy_equal_marginals)) legacy = legacy_equal_marginals
      res = homogeneity(tab, .false., legacy)
   end function

   pure function stuart_maxwell_ratings(ratings, legacy_equal_marginals) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Subject-by-two-rater categorical scores; NaN rows are omitted.
      logical, intent(in), optional :: legacy_equal_marginals !! Apply the R equal-marginal removal rule, default true.
      type(irr_result) :: res, base
      real(dp), allocatable :: tab(:, :)
      call contingency(ratings, tab, base)
      res = base
      if (base%status /= irr_ok) return
      res = stuart_maxwell_mh(tab, legacy_equal_marginals)
   end function

   pure function homogeneity(tab, bhapkar_correction, legacy) result(res)
      real(dp), intent(in) :: tab(:, :) !! Nonnegative finite square contingency counts.
      logical, intent(in) :: bhapkar_correction !! Subtract d*d'/n from the covariance for Bhapkar's test.
      logical, intent(in) :: legacy !! Remove all equal-marginal categories before Stuart-Maxwell computation.
      type(irr_result) :: res
      real(dp), allocatable :: x(:, :), rows(:), cols(:), d(:), cov(:, :), solution(:)
      integer, allocatable :: retained(:)
      integer :: n, i, j, info
      call init_result(res)
      if (.not. valid_table(tab)) return
      res%raters = 2
      rows = sum(tab, dim = 2)
      cols = sum(tab, dim = 1)
      n = size(tab, 1)
      if (legacy) then
         retained = pack([(i, i = 1,n)], rows /= cols)
      else
         ! Truly isolated perfectly agreeing categories carry no covariance information.
         retained = pack([(i, i = 1,n)], rows + cols - [(2.0_dp * tab(i, i), i = 1,n)] > 0.0_dp)
      end if
      if (bhapkar_correction) retained = [(i, i = 1,n)]
      x = tab(retained, retained)
      res%subjects = nint(sum(x))
      res%status = irr_degenerate
      n = size(x, 1) - 1
      if (n < 1) return
      rows = sum(x, dim = 2)
      cols = sum(x, dim = 1)
      d = rows(1:n) - cols(1:n)
      allocate(cov(n, n), solution(n))
      do j = 1, n
         do i = 1, n
            cov(i, j) = -x(i, j) - x(j, i)
            if (i == j) cov(i, j) = cov(i, j) + rows(i) + cols(i)
            if (bhapkar_correction) cov(i, j) = cov(i, j) - d(i) * d(j) / sum(x)
         end do
      end do
      call solve_system(cov, d, solution, info)
      if (info /= 0) return
      res%value = dot_product(d, solution)
      if (res%value < -100.0_dp * epsilon(1.0_dp)) return
      res%value = max(0.0_dp, res%value)
      res%statistic = res%value
      res%df1 = real(n, dp)
      res%p_value = r_pchisq(res%value, res%df1, lower_tail = .false.)
      res%status = irr_ok
   end function

   pure function rater_bias(tab) result(res)
      real(dp), intent(in) :: tab(:, :) !! Square nonnegative contingency counts with common ordered category labels.
      type(irr_result) :: res
      real(dp) :: upper, lower
      integer :: i, j
      call init_result(res)
      if (.not. valid_table(tab)) return
      res%subjects = nint(sum(tab))
      res%raters = 2
      upper = 0.0_dp
      lower = 0.0_dp
      do j = 1, size(tab, 2)
         do i = 1, j - 1
            upper = upper + tab(i, j)
            lower = lower + tab(j, i)
         end do
      end do
      res%status = irr_degenerate
      if (upper + lower == 0.0_dp) return
      res%value = upper / (upper + lower)
      res%statistic = (upper - lower) ** 2 / (upper + lower)
      res%df1 = 1.0_dp
      res%p_value = r_pchisq(res%statistic, 1.0_dp, lower_tail = .false.)
      res%status = irr_ok
   end function

   pure function rater_bias_ratings(ratings) result(res)
      real(dp), intent(in) :: ratings(:, :) !! Two rater score columns, explicitly distinguished from a contingency table.
      type(irr_result) :: res, base
      real(dp), allocatable :: tab(:, :)
      call contingency(ratings, tab, base)
      res = base
      if (base%status /= irr_ok) return
      res = rater_bias(tab)
      res%levels = base%levels
   end function
end module
