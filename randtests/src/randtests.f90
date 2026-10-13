! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from randtests 1.0.2 by Frederico Caeiro and Ayana Mateus.
module randtests
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only: dp
   use r_distributions, only: r_pnorm
   use r_special, only: r_regularized_beta
   use r_ordering, only: r_order
   use r_quantiles, only: r_median
   use randtests_distributions
   use randtests_permutations
   implicit none
   private
   integer, parameter, public :: rt_success = 0, rt_invalid = 1, rt_insufficient = 2
   integer, parameter, public :: rt_degenerate = 3, rt_limit = 4
   type, public :: randomness_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 0.0_dp
      real(dp) :: observed = 0.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: rvn = 0.0_dp
      integer :: n = 0
      integer :: n1 = 0
      integer :: n2 = 0
      integer :: status = rt_invalid
      character(len=128) :: message = ''
   end type randomness_result
   public :: dp, difference_sign_test, rank_test, turning_point_test
   public :: cox_stuart_test, runs_test, bartels_rank_test
   public :: druns, pruns, qruns, rruns, dbartelsrank, pbartelsrank, bartels_histogram
   public :: max_runs_total, max_exact_bartels_n, max_permutation_entries
   public :: permut, permut_apply, permutation_callback, randtests_aux
contains
   pure function difference_sign_test(x, alternative) result(test)
      real(dp), intent(in) :: x(:) !! Ordered observations; NaNs are omitted and infinities rejected.
      character(len=*), intent(in), optional :: alternative !! "two.sided" (default), "left.sided", "right.sided", or t/l/r.
      type(randomness_result) :: test
      real(dp), allocatable :: y(:)
      integer :: side, n, i, changes, rises
      call prepare_test(x, alternative, 2, y, side, test)
      if (test % status /= rt_success) return
      n = size(y)
      changes = 0
      rises = 0
      do i = 2, n
         if (y(i) /= y(i - 1)) changes = changes + 1
         if (y(i) > y(i - 1)) rises = rises + 1
      end do
      test % n = changes + 1
      test % observed = real(rises, dp)
      test % mean = real(changes, dp) / 2.0_dp
      test % variance = real(changes + 2, dp) / 12.0_dp
      call finish_normal(test, side)
   end function difference_sign_test

   pure function rank_test(x, alternative) result(test)
      real(dp), intent(in) :: x(:) !! Ordered observations; NaNs are omitted, infinities rejected, ties are not adjusted.
      character(len=*), intent(in), optional :: alternative !! "two.sided" (default), "left.sided", "right.sided", or t/l/r.
      type(randomness_result) :: test
      real(dp), allocatable :: y(:)
      integer :: side, i
      real(dp) :: n
      call prepare_test(x, alternative, 2, y, side, test)
      if (test % status /= rt_success) return
      test % observed = 0.0_dp
      do i = 1, size(y) - 1
         test % observed = test % observed + real(count(y(i + 1:) > y(i)), dp)
      end do
      n = real(size(y), dp)
      test % mean = n * (n - 1.0_dp) / 4.0_dp
      test % variance = n * (n - 1.0_dp) * (2.0_dp * n + 5.0_dp) / 72.0_dp
      call finish_normal(test, side)
   end function rank_test

   pure function turning_point_test(x, alternative) result(test)
      real(dp), intent(in) :: x(:) !! Ordered data; NaNs omitted, consecutive duplicates collapsed, infinities rejected.
      character(len=*), intent(in), optional :: alternative !! "two.sided" (default), "left.sided", "right.sided", or t/l/r.
      type(randomness_result) :: test
      real(dp), allocatable :: y(:), distinct(:)
      integer :: side, i, k
      real(dp) :: n
      call prepare_test(x, alternative, 3, y, side, test)
      if (test % status /= rt_success) return
      allocate (distinct(size(y)))
      k = 1
      distinct(1) = y(1)
      do i = 2, size(y)
         if (y(i) == distinct(k)) cycle
         k = k + 1
         distinct(k) = y(i)
      end do
      if (k < 3) then
         call fail_result(test, rt_insufficient, 'Need at least three observations after collapsing consecutive ties.')
         return
      end if
      test % n = k
      test % observed = 0.0_dp
      do i = 2, k - 1
         if ((distinct(i) > distinct(i - 1)) .neqv. (distinct(i + 1) > distinct(i))) &
            test % observed = test % observed + 1.0_dp
      end do
      n = real(k, dp)
      test % mean = 2.0_dp * (n - 2.0_dp) / 3.0_dp
      test % variance = (16.0_dp * n - 29.0_dp) / 90.0_dp
      call finish_normal(test, side)
   end function turning_point_test

   pure function cox_stuart_test(x, alternative) result(test)
      real(dp), intent(in) :: x(:) !! Ordered data; NaNs omitted, infinities rejected, odd middle observation dropped.
      character(len=*), intent(in), optional :: alternative !! "two.sided" (default), "left.sided", "right.sided", or t/l/r.
      type(randomness_result) :: test
      real(dp), allocatable :: y(:)
      integer :: side, half, offset, i, rises, pairs
      real(dp) :: lower, upper
      call prepare_test(x, alternative, 2, y, side, test)
      if (test % status /= rt_success) return
      half = size(y) / 2
      offset = size(y) - half
      pairs = 0
      rises = 0
      do i = 1, half
         if (y(i) /= y(i + offset)) pairs = pairs + 1
         if (y(i) < y(i + offset)) rises = rises + 1
      end do
      test % n = pairs
      test % observed = real(rises, dp)
      test % statistic = test % observed
      test % mean = real(pairs, dp) / 2.0_dp
      test % variance = real(pairs, dp) / 4.0_dp
      lower = binomial_half_cdf(rises, pairs)
      upper = binomial_half_cdf(pairs - rises, pairs)
      test % p_value = tail_probability(lower, upper, side)
   end function cox_stuart_test

   pure function runs_test(x, alternative, threshold, pvalue) result(test)
      real(dp), intent(in) :: x(:) !! Ordered observations; NaNs omitted, infinities rejected, threshold ties removed.
      character(len=*), intent(in), optional :: alternative !! "two.sided" (default), "left.sided", "right.sided", or t/l/r.
      real(dp), intent(in), optional :: threshold !! Finite classification level; default median of non-NaN observations.
      character(len=*), intent(in), optional :: pvalue !! P-value method: "normal" (default) or "exact".
      type(randomness_result) :: test
      real(dp), allocatable :: y(:), retained(:)
      real(dp) :: level, n, a, b, lower, upper
      integer :: side, i, runs
      logical :: exact
      call prepare_test(x, alternative, 2, y, side, test)
      if (test % status /= rt_success) return
      exact = .false.
      if (present(pvalue)) then
         if (pvalue /= 'normal' .and. pvalue /= 'exact') then
            call fail_result(test, rt_invalid, 'P-value method must be normal or exact.')
            return
         end if
         exact = pvalue == 'exact'
      end if
      level = r_median(y)
      if (present(threshold)) level = threshold
      if (.not. ieee_is_finite(level)) then
         call fail_result(test, rt_invalid, 'Threshold must be finite.')
         return
      end if
      retained = pack(y, y /= level)
      test % n1 = count(retained > level)
      test % n2 = count(retained < level)
      test % n = size(retained)
      if (test % n1 < 1 .or. test % n2 < 1) then
         call fail_result(test, rt_degenerate, 'Runs test needs observations on both sides of the threshold.')
         return
      end if
      if (test % n > max_runs_total) then
         call fail_result(test, rt_limit, 'Runs test supports at most 100000 retained observations.')
         return
      end if
      runs = 1
      do i = 2, size(retained)
         if ((retained(i) > level) .neqv. (retained(i - 1) > level)) runs = runs + 1
      end do
      n = real(test % n, dp)
      a = real(test % n1, dp)
      b = real(test % n2, dp)
      test % observed = real(runs, dp)
      test % mean = 1.0_dp + 2.0_dp * a * b / n
      test % variance = 2.0_dp * a * b * (2.0_dp * a * b - n) / (n ** 2 * (n - 1.0_dp))
      if (test % variance <= 0.0_dp) then
         if (exact) then
            test % statistic = 0.0_dp
            test % p_value = 1.0_dp
         else
            call fail_result(test, rt_degenerate, 'Normal runs approximation requires positive variance; use exact.')
         end if
         return
      end if
      call finish_normal(test, side)
      if (exact) then
         lower = pruns(test % observed, test % n1, test % n2)
         upper = pruns(test % observed - 1.0_dp, test % n1, test % n2, .false.)
         test % p_value = tail_probability(lower, upper, side)
      end if
   end function runs_test

   pure function bartels_rank_test(x, alternative, pvalue) result(test)
      real(dp), intent(in) :: x(:) !! Ordered data; NaNs omitted, infinities rejected, ties receive average ranks.
      character(len=*), intent(in), optional :: alternative !! "two.sided" (default), "left.sided", "right.sided", or t/l/r.
      character(len=*), intent(in), optional :: pvalue !! "normal" (default), "beta", "exact" (untied n<=9), or "auto".
      type(randomness_result) :: test
      real(dp), allocatable :: y(:), ranks(:)
      integer, allocatable :: order(:), histogram(:)
      integer :: side, i, j, n, status, score, total
      real(dp) :: nr, denominator, shape, lower, upper
      character(len=6) :: method
      logical :: tied
      call prepare_test(x, alternative, 3, y, side, test)
      if (test % status /= rt_success) return
      method = 'normal'
      if (present(pvalue)) then
         select case (pvalue)
         case ('normal', 'beta', 'exact', 'auto')
            method = pvalue
         case default
            call fail_result(test, rt_invalid, 'P-value method must be normal, beta, exact, or auto.')
            return
         end select
      end if
      n = size(y)
      if (method == 'auto') then
         method = 'normal'
         if (n <= 100) method = 'beta'
         if (n <= max_exact_bartels_n) method = 'exact'
      end if
      call r_order(y, order)
      allocate (ranks(n))
      tied = .false.
      i = 1
      do while (i <= n)
         j = i
         do while (j < n)
            if (y(order(j + 1)) /= y(order(i))) exit
            j = j + 1
         end do
         ranks(order(i:j)) = (real(i, dp) + real(j, dp)) / 2.0_dp
         tied = tied .or. j > i
         i = j + 1
      end do
      if (method == 'exact' .and. (n > max_exact_bartels_n .or. tied)) then
         call fail_result(test, rt_limit, 'Exact Bartels requires untied ranks and sample size at most nine.')
         return
      end if
      nr = real(n, dp)
      denominator = sum((ranks - (nr + 1.0_dp) / 2.0_dp) ** 2)
      if (denominator <= 0.0_dp) then
         call fail_result(test, rt_degenerate, 'Bartels test is undefined when all observations are equal.')
         return
      end if
      test % observed = sum((ranks(2:) - ranks(:n - 1)) ** 2)
      test % rvn = test % observed / denominator
      test % mean = 2.0_dp
      test % variance = 4.0_dp * (nr - 2.0_dp) * (5.0_dp * nr ** 2 - 2.0_dp * nr - 9.0_dp) / &
                        (5.0_dp * nr * (nr + 1.0_dp) * (nr - 1.0_dp) ** 2)
      test % statistic = (test % rvn - test % mean) / sqrt(test % variance)
      select case (method)
      case ('normal')
         lower = r_pnorm(test % statistic)
         upper = r_pnorm(test % statistic, lower_tail=.false.)
      case ('beta')
         shape = 2.0_dp / test % variance - 0.5_dp
         lower = r_regularized_beta(test % rvn / 4.0_dp, shape, shape)
         upper = r_regularized_beta(1.0_dp - test % rvn / 4.0_dp, shape, shape)
      case ('exact')
         call bartels_histogram(n, histogram, status)
         score = nint(test % observed)
         total = sum(histogram)
         lower = real(sum(histogram(:score)), dp) / real(total, dp)
         upper = real(sum(histogram(score:)), dp) / real(total, dp)
      end select
      test % p_value = tail_probability(lower, upper, side)
   end function bartels_rank_test

   pure subroutine prepare_test(x, alternative, minimum_n, y, side, test)
      real(dp), intent(in) :: x(:) !! Raw ordered observations, permitting NaNs but not infinities.
      character(len=*), intent(in), optional :: alternative !! Tail selection; absent means two-sided.
      integer, intent(in) :: minimum_n !! Minimum finite sample size for this test before further preprocessing.
      real(dp), allocatable, intent(out) :: y(:) !! Input observations with NaNs removed in original order.
      integer, intent(out) :: side !! Tail selector: zero two-sided, -1 lower, +1 upper.
      type(randomness_result), intent(out) :: test !! Initialized result, with status and message for invalid input.
      call fail_result(test, rt_invalid, 'Invalid alternative; use two.sided, left.sided, or right.sided.')
      side = 0
      if (present(alternative)) then
         select case (alternative)
         case ('two.sided', 't')
            side = 0
         case ('left.sided', 'l')
            side = -1
         case ('right.sided', 'r')
            side = 1
         case default
            return
         end select
      end if
      y = pack(x,.not. ieee_is_nan(x))
      if (any(.not. ieee_is_finite(y))) then
         call fail_result(test, rt_invalid, 'Infinite observations are not supported.')
         return
      end if
      if (size(y) < minimum_n) then
         call fail_result(test, rt_insufficient, 'Too few observations after omitting NaNs.')
         return
      end if
      test = randomness_result()
      test % status = rt_success
      test % n = size(y)
   end subroutine prepare_test

   pure subroutine fail_result(test, status, message)
      type(randomness_result), intent(out) :: test !! Error result whose real-valued fields are all quiet NaNs.
      integer, intent(in) :: status !! Nonzero error classification from the public rt_* constants.
      character(len=*), intent(in) :: message !! Human-readable reason for the rejected input or unsupported operation.
      test = randomness_result()
      test % status = status
      test % message = message
      test % statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      test % p_value = test % statistic
      test % observed = test % statistic
      test % mean = test % statistic
      test % variance = test % statistic
      test % rvn = test % statistic
   end subroutine fail_result

   pure subroutine finish_normal(test, side)
      type(randomness_result), intent(inout) :: test !! Result with observed, mean, and positive variance already set.
      integer, intent(in) :: side !! Tail selector: zero two-sided, -1 lower, +1 upper.
      test % statistic = (test % observed - test % mean) / sqrt(test % variance)
      test % p_value = tail_probability(r_pnorm(test % statistic), &
                                        r_pnorm(test % statistic, lower_tail=.false.), side)
   end subroutine finish_normal

   pure elemental real(dp) function tail_probability(lower, upper, side) result(p)
      real(dp), intent(in) :: lower !! Lower tail, including the observed statistic for discrete distributions.
      real(dp), intent(in) :: upper !! Upper tail, including the observed statistic for discrete distributions.
      integer, intent(in) :: side !! Zero gives doubled smaller tail capped at one; -1 lower, +1 upper.
      if (side == -1) then
         p = lower
      else if (side == 1) then
         p = upper
      else
         p = min(1.0_dp, 2.0_dp * min(lower, upper))
      end if
   end function tail_probability

   pure elemental real(dp) function binomial_half_cdf(k, n) result(p)
      integer, intent(in) :: k !! Integer upper cutoff for a binomial(n,0.5) lower tail.
      integer, intent(in) :: n !! Nonnegative number of Bernoulli trials.
      if (k < 0) then
         p = 0.0_dp
      else if (k >= n) then
         p = 1.0_dp
      else
         p = r_regularized_beta(0.5_dp, real(n - k, dp), real(k + 1, dp))
      end if
   end function binomial_half_cdf
end module randtests
