! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from exactRankTests by Torsten Hothorn and Kurt Hornik.
module exact_rank_tests
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_negative_inf
   use r_kinds, only : dp
   use r_sorting, only : r_average_ranks
   use r_ordering, only : r_order, r_sort_values_in_place
   use r_distributions, only : r_pnorm
   use exact_rank_distribution
   use exact_rank_scores
   implicit none
   private
   type, public :: rank_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: point_probability = 0.0_dp
      real(dp) :: expected = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: z = 0.0_dp
      real(dp) :: estimate = 0.0_dp
      real(dp) :: conf_int(2) = 0.0_dp
      logical :: exact = .false.
      logical :: has_conf_int = .false.
      integer :: status = rank_invalid
   end type rank_test_result
   public :: dp, rank_ok, rank_invalid, rank_limit, rank_allocation, rank_numeric
   public :: wilcox_exact, ansari_exact, perm_test
   public :: rank_distribution, permutation_distribution, simulate_distribution, dperm, pperm, qperm, rperm
   public :: cscores, cscores_factor, cscores_surv, irank
contains
   pure subroutine initialize(result)
      type(rank_test_result), intent(out) :: result !! Reset test record to an invalid, IEEE-NaN result.
      result%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      result%p_value = result%statistic
      result%point_probability = result%statistic
      result%expected = result%statistic
      result%variance = result%statistic
      result%z = result%statistic
      result%estimate = result%statistic
      result%conf_int = result%statistic
   end subroutine initialize

   pure subroutine score_test(scores, m, statistic, paired, exact, result, alternative, tolerance)
      real(dp), intent(in) :: scores(:) !! Finite rank or data scores defining the permutation population.
      integer, intent(in) :: m !! Size of the selected first group, or full size for sign flips.
      real(dp), intent(in) :: statistic !! Observed uncentered sum of selected scores.
      logical, intent(in) :: paired !! True for sign-flip inference, false for fixed-size sampling.
      logical, intent(in) :: exact !! True requests the conditional lattice distribution; false requests normal approximation.
      type(rank_test_result), intent(out) :: result !! Numerical statistic, moments, probabilities, and status.
      character(*), intent(in), optional :: alternative !! 'two.sided' (default), 'less', or 'greater'.
      real(dp), intent(in), optional :: tolerance !! Lattice approximation tolerance, default 0.01.
      type(rank_distribution) :: distribution
      character(9) :: alt
      real(dp) :: n, lower, upper
      call initialize(result)
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) return
         alt = alternative
      end if
      if (alt /= 'two.sided' .and. alt /= 'less' .and. alt /= 'greater') return
      if (m < 1 .or. m > size(scores)) return
      if (any(.not. ieee_is_finite(scores)) .or. .not. ieee_is_finite(statistic)) return
      n = real(size(scores), dp)
      result%statistic = statistic
      result%exact = exact
      if (paired) then
         result%expected = sum(0.5_dp * scores)
         result%variance = sum((0.5_dp * scores) ** 2)
      else
         result%expected = sum(scores / n) * real(m, dp)
         result%variance = 0.0_dp
         if (size(scores) > 1) result%variance = &
            (real(m, dp) / n) * ((n - real(m, dp)) / (n - 1.0_dp)) * &
            sum((scores - sum(scores / n)) ** 2)
      end if
      if (.not. ieee_is_finite(result%variance)) then
         result%status = rank_numeric
         return
      end if
      if (result%variance > 0.0_dp) result%z = (statistic - result%expected) / sqrt(result%variance)
      if (exact) then
         call permutation_distribution(scores, m, distribution, result%status, paired, tolerance = tolerance)
         if (result%status /= rank_ok) return
         result%p_value = pperm(statistic, distribution, alt)
         result%point_probability = dperm(statistic, distribution)
      else
         if (result%variance <= 0.0_dp) then
            result%status = rank_numeric
            return
         end if
         lower = r_pnorm(result%z)
         upper = r_pnorm(result%z, lower_tail = .false.)
         select case (alt)
         case ('two.sided')
            result%p_value = 2.0_dp * min(lower, upper)
         case ('less')
            result%p_value = lower
         case ('greater')
            result%p_value = upper
         end select
         result%status = rank_ok
      end if
   end subroutine score_test

   pure subroutine wilcox_exact(x, result, y, paired, alternative, mu, exact, conf_level)
      real(dp), intent(in) :: x(:) !! First sample or one-sample observations; nonfinite observations are omitted.
      type(rank_test_result), intent(out) :: result !! V for signed ranks or W for rank sums, p-value, moments, and optional interval.
      real(dp), intent(in), optional :: y(:) !! Second sample; paired mode requires the same length as x.
      logical, intent(in), optional :: paired !! Paired differences x-y when y is present; default false.
      character(*), intent(in), optional :: alternative !! 'two.sided' (default), 'less', or 'greater'.
      real(dp), intent(in), optional :: mu !! Finite null location shift, default zero.
      logical, intent(in), optional :: exact !! Default true when each effective sample has fewer than 50 observations.
      real(dp), intent(in), optional :: conf_level !! Request an exact shift interval and estimate; must be in (0,1) and exact=true.
      real(dp), allocatable :: first(:), second(:), original(:), ranks(:), values(:)
      real(dp) :: shift, observed, offset
      integer :: m, n
      logical :: paired_data, sign_flip, use_exact
      call initialize(result)
      shift = 0.0_dp
      if (present(mu)) shift = mu
      if (.not. ieee_is_finite(shift)) return
      paired_data = .false.
      if (present(paired)) paired_data = paired
      sign_flip = .not. present(y) .or. paired_data
      if (paired_data .and. .not. present(y)) return
      if (paired_data) then
         if (size(x) /= size(y)) return
         first = pack(x, ieee_is_finite(x) .and. ieee_is_finite(y)) - &
            pack(y, ieee_is_finite(x) .and. ieee_is_finite(y))
      else
         first = pack(x, ieee_is_finite(x))
      end if
      if (size(first) < 1) return
      if (any(.not. ieee_is_finite(first - shift))) return
      original = first
      first = first - shift
      offset = 0.0_dp
      if (sign_flip) then
         first = pack(first, first /= 0.0_dp)
         m = size(first)
         if (m < 1) return
         call r_average_ranks(abs(first), ranks)
         observed = sum(ranks, mask = first > 0.0_dp)
         use_exact = m < 50
      else
         second = pack(y, ieee_is_finite(y))
         m = size(first)
         n = size(second)
         if (n < 1) return
         values = [first, second]
         call r_average_ranks(values, ranks)
         observed = sum(ranks(:m))
         offset = real(m, dp) * real(m + 1, dp) / 2.0_dp
         use_exact = m < 50 .and. n < 50
      end if
      if (present(exact)) use_exact = exact
      call score_test(ranks, m, observed, sign_flip, use_exact, result, alternative)
      result%statistic = result%statistic - offset
      result%expected = result%expected - offset
      if (result%status /= rank_ok) return
      if (present(conf_level)) then
         if (.not. ieee_is_finite(conf_level) .or. conf_level <= 0.0_dp .or. conf_level >= 1.0_dp) then
            result%status = rank_invalid
            return
         end if
         if (.not. use_exact) then
            result%status = rank_invalid
            return
         end if
         first = first + shift
         if (sign_flip) then
            call wilcox_interval(first, ranks, result, conf_level, alternative)
         else
            call wilcox_interval(original, ranks, result, conf_level, alternative, second)
         end if
      end if
   end subroutine wilcox_exact

   pure subroutine wilcox_interval(x, tested_ranks, result, level, alternative, y)
      real(dp), intent(in) :: x(:) !! Retained first sample before null-shift subtraction.
      real(dp), intent(in) :: tested_ranks(:) !! Ranks used for the original test, including its null shift.
      type(rank_test_result), intent(inout) :: result !! Successful exact test record to receive confidence interval and estimate.
      real(dp), intent(in) :: level !! Confidence level strictly between zero and one.
      character(*), intent(in), optional :: alternative !! Alternative controlling one- or two-sided interval endpoints.
      real(dp), intent(in), optional :: y(:) !! Retained second sample; absence means signed-rank inference.
      real(dp), allocatable :: diffs(:), w(:), ranks(:), values(:), base_ranks(:), shifted(:)
      type(rank_distribution) :: distribution
      integer :: i, j, k, m, n, length
      real(dp) :: alpha, lower_quantile, upper_quantile, offset, mean_rank, low, high
      logical :: ties, sign_flip
      character(9) :: alt
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) return
         alt = alternative
      end if
      m = size(x)
      sign_flip = .not. present(y)
      offset = 0.0_dp
      if (sign_flip) then
         if (real(m, dp) * real(m + 1, dp) / 2.0_dp > 1000000.0_dp) then
            result%status = rank_limit
            return
         end if
         length = m * (m + 1) / 2
         allocate(diffs(length))
         k = 0
         do i = 1, m
            do j = i, m
               k = k + 1
               diffs(k) = x(i) / 2.0_dp + x(j) / 2.0_dp
            end do
         end do
         base_ranks = tested_ranks
      else
         n = size(y)
         if (real(m, dp) * real(n, dp) > 1000000.0_dp) then
            result%status = rank_limit
            return
         end if
         length = m * n
         allocate(diffs(length))
         k = 0
         do i = 1, m
            do j = 1, n
               k = k + 1
               diffs(k) = x(i) - y(j)
            end do
         end do
         values = [x, y]
         call r_average_ranks(values, base_ranks)
         offset = real(m, dp) * real(m + 1, dp) / 2.0_dp
      end if
      if (any(.not. ieee_is_finite(diffs))) then
         result%status = rank_numeric
         return
      end if
      call r_sort_values_in_place(diffs)
      values = tested_ranks
      call r_sort_values_in_place(values)
      ties = any(values(2:) == values(:size(values) - 1))
      allocate(w(length))
      if (ties) then
         do k = 1, length
            shifted = x - diffs(k)
            if (sign_flip) then
               call r_average_ranks(abs(shifted), ranks)
               w(k) = sum(ranks, mask = shifted > 0.0_dp)
            else
               values = [shifted, y]
               call r_average_ranks(values, ranks)
               w(k) = sum(ranks(:m)) - offset
            end if
         end do
      else
         w = [(real(length - k + 1, dp), k = 1, length)]
      end if
      call permutation_distribution(base_ranks, m, distribution, result%status, paired = sign_flip)
      if (result%status /= rank_ok) return
      alpha = 1.0_dp - level
      if (alt == 'two.sided') alpha = alpha / 2.0_dp
      lower_quantile = qperm(alpha, distribution) - offset
      upper_quantile = qperm(1.0_dp - alpha, distribution) - offset
      low = minval(diffs)
      high = maxval(diffs)
      if (upper_quantile < maxval(w)) low = maxval(diffs, mask = w > upper_quantile)
      if (lower_quantile > minval(w)) high = minval(diffs, mask = w <= lower_quantile)
      if (alt == 'less') low = ieee_value(0.0_dp, ieee_negative_inf)
      if (alt == 'greater') high = ieee_value(0.0_dp, ieee_positive_inf)
      result%conf_int = [low, high]
      mean_rank = distribution%expectation - offset
      if (any(w <= real(ceiling(mean_rank), dp)) .and. any(w > mean_rank)) then
         result%estimate = minval(diffs, mask = w <= real(ceiling(mean_rank), dp)) / 2.0_dp + &
            maxval(diffs, mask = w > mean_rank) / 2.0_dp
      end if
      result%has_conf_int = .true.
   end subroutine wilcox_interval

   pure subroutine ansari_exact(x, y, result, alternative, exact, conf_level)
      real(dp), intent(in) :: x(:) !! First finite sample; NaNs are omitted.
      real(dp), intent(in) :: y(:) !! Second finite sample; NaNs are omitted.
      type(rank_test_result), intent(out) :: result !! Ansari-Bradley statistic and exact or asymptotic p-value.
      character(*), intent(in), optional :: alternative !! 'two.sided' (default), 'less', or 'greater', following upstream conventions.
      logical, intent(in), optional :: exact !! Default true when both effective sample sizes are below 50.
      real(dp), intent(in), optional :: conf_level !! Exact scale-ratio interval level in (0,1); finite nonzero samples required.
      real(dp), allocatable :: first(:), second(:), ranks(:), values(:)
      integer :: m, n, status
      logical :: use_exact
      character(9) :: alt
      call initialize(result)
      first = pack(x, .not. ieee_is_nan(x))
      second = pack(y, .not. ieee_is_nan(y))
      m = size(first)
      n = size(second)
      if (m < 1 .or. n < 1) return
      values = [first, second]
      call cscores(values, ranks, status, 'AnsariBradley')
      if (status /= rank_ok) return
      use_exact = m < 50 .and. n < 50
      if (present(exact)) use_exact = exact
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) return
         alt = alternative
      end if
      if (.not. use_exact) then
         ! exactRankTests reverses its asymptotic one-sided Ansari convention.
         if (alt == 'less') then
            alt = 'greater'
         else if (alt == 'greater') then
            alt = 'less'
         end if
      end if
      call score_test(ranks, m, sum(ranks(:m)), .false., use_exact, result, alt)
      if (result%status /= rank_ok) return
      if (present(conf_level)) then
         if (.not. ieee_is_finite(conf_level) .or. conf_level <= 0.0_dp .or. conf_level >= 1.0_dp) then
            result%status = rank_invalid
            return
         end if
         if (.not. use_exact .or. any(first == 0.0_dp) .or. any(second == 0.0_dp)) then
            result%status = rank_invalid
            return
         end if
         call ansari_interval(first, second, ranks, result, conf_level, alt)
      end if
   end subroutine ansari_exact

   pure subroutine ansari_interval(x, y, scores, result, level, alternative)
      real(dp), intent(in) :: x(:) !! Finite nonzero first-sample observations.
      real(dp), intent(in) :: y(:) !! Finite nonzero second-sample observations.
      real(dp), intent(in) :: scores(:) !! Ansari-Bradley scores of the unscaled pooled observations.
      type(rank_test_result), intent(inout) :: result !! Test result to receive the exact scale interval and estimate.
      real(dp), intent(in) :: level !! Confidence level strictly between zero and one.
      character(*), intent(in) :: alternative !! 'two.sided', 'less', or 'greater'.
      type(rank_distribution) :: distribution
      real(dp), allocatable :: ratio(:), sigma(:), steps(:), values(:), ranks(:), u(:), l(:)
      real(dp) :: alpha, evaluation, qu, ql, lower, upper, best, mid
      integer :: i, j, k, m, n, count_ratio
      m = size(x)
      n = size(y)
      if (real(m, dp) * real(n, dp) > 1000000.0_dp) then
         result%status = rank_limit
         return
      end if
      allocate(ratio(m * n))
      k = 0
      do i = 1, m
         do j = 1, n
            k = k + 1
            if (abs(y(j)) < 1.0_dp) then
               if (abs(x(i)) > huge(1.0_dp) * abs(y(j))) then
                  result%status = rank_numeric
                  return
               end if
            end if
            ratio(k) = x(i) / y(j)
         end do
      end do
      sigma = pack(ratio, ratio > 0.0_dp)
      count_ratio = size(sigma)
      result%has_conf_int = .true.
      if (count_ratio == 0) return
      call r_sort_values_in_place(sigma)
      allocate(steps(count_ratio))
      do k = 1, count_ratio
         if (k < count_ratio) then
            evaluation = sigma(k) / 2.0_dp + sigma(k + 1) / 2.0_dp
         else
            if (sigma(k) > huge(1.0_dp) / 2.01_dp) then
               result%status = rank_numeric
               return
            end if
            evaluation = sigma(k) * 2.01_dp
         end if
         if (evaluation <= 0.0_dp) then
            result%status = rank_numeric
            return
         end if
         if (evaluation < 1.0_dp) then
            if (maxval(abs(x)) > huge(1.0_dp) * evaluation) then
               result%status = rank_numeric
               return
            end if
         end if
         values = [x / evaluation, y]
         if (any(.not. ieee_is_finite(values))) then
            result%status = rank_numeric
            return
         end if
         call r_average_ranks(values, ranks)
         ranks = min(ranks, real(m + n + 1, dp) - ranks)
         steps(k) = sum(ranks(:m))
      end do
      call permutation_distribution(scores, m, distribution, result%status, paired = .false.)
      if (result%status /= rank_ok) return
      alpha = (1.0_dp - level) / 2.0_dp
      if (alternative /= 'two.sided') alpha = 1.0_dp - level
      qu = qperm(alpha, distribution)
      ql = qperm(1.0_dp - alpha, distribution)
      u = steps - qu
      l = steps - ql
      if (any(u >= 0.0_dp) .and. any(l > 0.0_dp)) then
         best = minval(u, mask = u >= 0.0_dp)
         lower = minval(sigma, mask = u == best)
         best = minval(l, mask = l > 0.0_dp)
         upper = maxval(sigma, mask = l == best)
         if (lower > upper) then
            u = steps - ql
            l = steps - qu
            if (any(u >= 0.0_dp) .and. any(l > 0.0_dp)) then
               best = minval(u, mask = u >= 0.0_dp)
               lower = minval(sigma, mask = u == best)
               best = minval(l, mask = l > 0.0_dp)
               upper = maxval(sigma, mask = l == best)
            else
               lower = ieee_value(0.0_dp, ieee_quiet_nan)
               upper = lower
            end if
         end if
         result%conf_int = [lower, upper]
      end if
      if (alternative == 'less') result%conf_int(1) = 0.0_dp
      if (alternative == 'greater') result%conf_int(2) = ieee_value(0.0_dp, ieee_positive_inf)
      mid = qperm(0.5_dp, distribution)
      if (any(steps <= mid) .and. any(steps > mid)) then
         result%estimate = maxval(sigma, mask = steps <= mid) / 2.0_dp + minval(sigma, mask = steps > mid) / 2.0_dp
      end if
   end subroutine ansari_interval

   pure subroutine perm_test(x, result, y, paired, alternative, mu, exact, tolerance)
      real(dp), intent(in) :: x(:) !! First sample or one-sample differences; NaNs are omitted.
      type(rank_test_result), intent(out) :: result !! Permutation sum statistic, moments, point mass, and p-value.
      real(dp), intent(in), optional :: y(:) !! Optional second sample, same length as x for paired inference.
      logical, intent(in), optional :: paired !! Use x-y differences; absence of y always selects sign flips.
      character(*), intent(in), optional :: alternative !! 'two.sided' (default), 'less', or 'greater'.
      real(dp), intent(in), optional :: mu !! Finite null location shift, default zero.
      logical, intent(in), optional :: exact !! Defaults to effective sample sizes <=50.
      real(dp), intent(in), optional :: tolerance !! Explicit lattice tolerance; omitting it enables upstream score remapping.
      real(dp), allocatable :: first(:), second(:), scores(:), mapped(:)
      integer, allocatable :: order(:)
      integer :: m, n, i, negative
      real(dp) :: shift, observed, value, a
      logical :: paired_data, sign_flip, use_exact
      call initialize(result)
      paired_data = .false.
      if (present(paired)) paired_data = paired
      sign_flip = paired_data .or. .not. present(y)
      shift = 0.0_dp
      if (present(mu)) shift = mu
      if (.not. ieee_is_finite(shift)) return
      if (paired_data .and. present(y)) then
         if (size(x) /= size(y)) return
         first = pack(x, .not. ieee_is_nan(x) .and. .not. ieee_is_nan(y)) - &
            pack(y, .not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
      else
         first = pack(x, .not. ieee_is_nan(x))
      end if
      first = first - shift
      if (any(.not. ieee_is_finite(first))) return
      m = size(first)
      if (m < 1) return
      use_exact = m <= 50
      if (.not. sign_flip) then
         second = pack(y, .not. ieee_is_nan(y))
         n = size(second)
         if (n < 1 .or. any(.not. ieee_is_finite(second))) return
         use_exact = use_exact .and. n <= 50
      end if
      if (present(exact)) use_exact = exact
      if (sign_flip) then
         if (use_exact .and. .not. present(tolerance)) then
            if (any(first /= aint(first))) then
               call r_order(first, order)
               negative = count(first <= 0.0_dp)
               allocate(mapped(m))
               do i = 1, m
                  mapped(order(i)) = real(i - 1 - negative, dp)
               end do
               first = mapped
            end if
         end if
         observed = sum(first, mask = first > 0.0_dp)
         scores = abs(first)
      else
         scores = [first, second]
         if (use_exact .and. .not. present(tolerance)) then
            if (any(scores /= aint(scores))) then
               scores = scores - minval(scores)
               if (maxval(scores) > 0.0_dp) scores = scores / maxval(scores) * real(m + n, dp)
               do i = 1, m + n
                  value = scores(i)
                  a = aint(value)
                  if (value - a > 0.5_dp .or. (value - a == 0.5_dp .and. modulo(a, 2.0_dp) /= 0.0_dp)) &
                     a = a + 1.0_dp
                  scores(i) = a
               end do
            end if
         end if
         observed = sum(scores(:m))
      end if
      call score_test(scores, m, observed, sign_flip, use_exact, result, alternative, tolerance)
   end subroutine perm_test
end module exact_rank_tests
