! SPDX-License-Identifier: GPL-2.0-only
! Derived from matrixTests 0.2.3.1 by Karolis Koncevicius.
module matrixtests_kernels
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   use matrixtests_types
   use r_descriptive, only : r_mean, r_variance
   use r_distributions, only : r_pt, r_qt, r_pf, r_qf, r_pnorm, r_qnorm, r_pchisq
   use r_sorting, only : r_average_ranks
   use matrixstats, only : row_medians
   use exact_rank_distribution, only : rank_distribution, permutation_distribution, pperm, rank_ok
   use nortest, only : normality_result, ad_test
   use moments, only : moments_test_result, jarque_test, skewness, kurtosis
   use la_lapack_lsq, only : la_dgelsy
   implicit none
   private
   public :: t_kernel, f_kernel, correlation_kernel, group_kernel, wilcoxon_kernel
   public :: ad_kernel, jarque_kernel, cosinor_kernel, valid_alternative
contains
   pure elemental logical function valid_alternative(alternative) result(valid)
      character(*), intent(in) :: alternative !! Exact alternative name: two.sided, less, or greater.
      valid = alternative == 'two.sided' .or. alternative == 'less' .or. alternative == 'greater'
   end function valid_alternative

   pure function t_kernel(x, y, mode, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:) !! First sample, with IEEE NaNs omitted.
      real(dp), intent(in) :: y(:) !! Second sample; paired modes require the same length as x.
      integer, intent(in) :: mode !! 1 one-sample, 2 equal-variance, 3 Welch, or 4 paired t test.
      real(dp), intent(in) :: null !! Finite null mean or null difference of means.
      character(*), intent(in) :: alternative !! two.sided, less, or greater.
      real(dp), intent(in) :: conf_level !! Confidence level in [0,1]; NaN disables confidence intervals.
      type(test_result) :: r
      real(dp), allocatable :: a(:), b(:), d(:)
      real(dp) :: se, df, estimate, scale, va, vb, q
      integer :: nx, ny, nd
      r = empty_result()
      r%null = null
      r%alternative = alternative
      r%conf_level = conf_level
      if (.not. valid_alternative(alternative) .or. .not. ieee_is_finite(null)) return
      if (.not. ieee_is_nan(conf_level)) then
         if (conf_level < 0.0_dp .or. conf_level > 1.0_dp) return
      end if
      a = pack(x, .not. ieee_is_nan(x))
      b = pack(y, .not. ieee_is_nan(y))
      nx = size(a)
      ny = size(b)
      r%obs = nx
      r%obs_x = nx
      r%obs_y = ny
      r%obs_tot = nx + ny
      r%mean_x = r_mean(a)
      r%mean_y = r_mean(b)
      r%var_x = r_variance(a)
      r%var_y = r_variance(b)
      r%mean = r%mean_x
      r%var = r%var_x
      r%mean_diff = r%mean_x - r%mean_y
      r%status = mt_insufficient
      if (mode == 1) then
         if (nx < 2) return
         estimate = r%mean
         se = sqrt(r%var / real(nx, dp))
         df = real(nx - 1, dp)
         scale = abs(estimate)
      else if (mode == 4) then
         if (size(x) /= size(y)) then
            r%status = mt_invalid
            return
         end if
         d = x - y
         d = pack(d, .not. ieee_is_nan(d))
         nd = size(d)
         r%obs_paired = nd
         r%mean_diff = r_mean(d)
         r%var_diff = r_variance(d)
         if (nd < 2) return
         estimate = r%mean_diff
         se = sqrt(r%var_diff / real(nd, dp))
         df = real(nd - 1, dp)
         scale = abs(estimate)
      else
         estimate = r%mean_diff
         scale = max(abs(r%mean_x), abs(r%mean_y))
         if (mode == 2) then
            if (nx < 1 .or. ny < 1 .or. nx + ny < 3) return
            df = real(nx + ny - 2, dp)
            r%var_pooled = 0.0_dp
            if (nx > 1) r%var_pooled = r%var_pooled + (real(nx - 1, dp) / df) * r%var_x
            if (ny > 1) r%var_pooled = r%var_pooled + (real(ny - 1, dp) / df) * r%var_y
            se = sqrt(r%var_pooled * (1.0_dp / real(nx, dp) + 1.0_dp / real(ny, dp)))
         else
            if (nx < 2 .or. ny < 2) return
            va = r%var_x / real(nx, dp)
            vb = r%var_y / real(ny, dp)
            se = sqrt(va + vb)
            df = 1.0_dp / ((va / (va + vb)) ** 2 / real(nx - 1, dp) + &
               (vb / (va + vb)) ** 2 / real(ny - 1, dp))
         end if
      end if
      r%status = mt_degenerate
      if (se <= 10.0_dp * epsilon(1.0_dp) * scale) return
      r%status = mt_nonfinite
      if (.not. ieee_is_finite(se) .or. .not. ieee_is_finite(estimate)) return
      if (.not. ieee_is_finite(df)) return
      r%stderr = se
      r%df = df
      r%statistic = (estimate - null) / se
      r%pvalue = student_tail(r%statistic, df, alternative)
      if (.not. ieee_is_nan(conf_level)) then
         if (alternative == 'two.sided') then
            q = r_qt((1.0_dp - conf_level) / 2.0_dp, df, lower_tail = .false.)
            r%conf_low = estimate - q * se
            r%conf_high = estimate + q * se
         else if (alternative == 'less') then
            r%conf_low = ieee_value(0.0_dp, ieee_negative_inf)
            r%conf_high = estimate + r_qt(conf_level, df) * se
         else
            r%conf_low = estimate - r_qt(conf_level, df) * se
            r%conf_high = ieee_value(0.0_dp, ieee_positive_inf)
         end if
      end if
      r%status = mt_ok
   end function t_kernel

   pure elemental real(dp) function student_tail(t, df, alternative) result(p)
      real(dp), intent(in) :: t !! Signed Student t statistic, possibly infinite.
      real(dp), intent(in) :: df !! Positive Student distribution degrees of freedom.
      character(*), intent(in) :: alternative !! Validated tail selection.
      if (alternative == 'two.sided') then
         p = 2.0_dp * r_pt(-abs(t), df)
      else
         p = r_pt(t, df, lower_tail = alternative == 'less')
      end if
   end function student_tail

   pure function f_kernel(x, y, null, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:) !! First sample; NaNs and infinities are removed as upstream.
      real(dp), intent(in) :: y(:) !! Second sample; NaNs and infinities are removed as upstream.
      real(dp), intent(in) :: null !! Positive finite null variance ratio.
      character(*), intent(in) :: alternative !! two.sided, less, or greater.
      real(dp), intent(in) :: conf_level !! Confidence level in [0,1]; NaN disables intervals.
      type(test_result) :: r
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: pl, pu, beta
      r = empty_result()
      r%null = null
      r%alternative = alternative
      r%conf_level = conf_level
      if (.not. valid_alternative(alternative) .or. .not. ieee_is_finite(null)) return
      if (null <= 0.0_dp) return
      if (.not. ieee_is_nan(conf_level)) then
         if (conf_level < 0.0_dp .or. conf_level > 1.0_dp) return
      end if
      a = pack(x, ieee_is_finite(x))
      b = pack(y, ieee_is_finite(y))
      r%obs_x = size(a)
      r%obs_y = size(b)
      r%obs_tot = size(a) + size(b)
      r%var_x = r_variance(a)
      r%var_y = r_variance(b)
      r%status = mt_insufficient
      if (size(a) < 2 .or. size(b) < 2) return
      r%status = mt_nonfinite
      if (.not. ieee_is_finite(r%var_x) .or. .not. ieee_is_finite(r%var_y)) return
      r%status = mt_degenerate
      if (r%var_x == 0.0_dp .and. r%var_y == 0.0_dp) return
      if (r%var_y == 0.0_dp) then
         r%var_ratio = ieee_value(0.0_dp, ieee_positive_inf)
      else
         r%var_ratio = r%var_x / r%var_y
      end if
      r%df_num = real(size(a) - 1, dp)
      r%df_denom = real(size(b) - 1, dp)
      r%statistic = r%var_ratio / null
      pl = r_pf(r%statistic, r%df_num, r%df_denom)
      pu = r_pf(r%statistic, r%df_num, r%df_denom, lower_tail = .false.)
      select case (alternative)
      case ('two.sided')
         r%pvalue = 2.0_dp * min(pl, pu)
      case ('less')
         r%pvalue = pl
      case ('greater')
         r%pvalue = pu
      end select
      if (.not. ieee_is_nan(conf_level)) then
         if (alternative == 'two.sided') then
            beta = (1.0_dp - conf_level) / 2.0_dp
            r%conf_low = r%var_ratio / r_qf(beta, r%df_num, r%df_denom, lower_tail = .false.)
            r%conf_high = r%var_ratio / r_qf(beta, r%df_num, r%df_denom)
         else if (alternative == 'less') then
            r%conf_low = 0.0_dp
            r%conf_high = r%var_ratio / r_qf(conf_level, r%df_num, r%df_denom, lower_tail = .false.)
         else
            r%conf_low = r%var_ratio / r_qf(conf_level, r%df_num, r%df_denom)
            r%conf_high = ieee_value(0.0_dp, ieee_positive_inf)
         end if
      end if
      r%status = mt_ok
   end function f_kernel

   pure function correlation_kernel(x, y, alternative, conf_level) result(r)
      real(dp), intent(in) :: x(:) !! First sample, paired with y; missing pairs are omitted.
      real(dp), intent(in) :: y(:) !! Second sample of the same length as x.
      character(*), intent(in) :: alternative !! two.sided, less, or greater, against zero correlation.
      real(dp), intent(in) :: conf_level !! Confidence level in [0,1]; NaN disables intervals.
      type(test_result) :: r
      real(dp), allocatable :: a(:), b(:)
      logical, allocatable :: keep(:)
      real(dp) :: sx, sy, z, q, se
      r = empty_result()
      r%null = 0.0_dp
      r%alternative = alternative
      r%conf_level = conf_level
      if (size(x) /= size(y) .or. .not. valid_alternative(alternative)) return
      if (.not. ieee_is_nan(conf_level)) then
         if (conf_level < 0.0_dp .or. conf_level > 1.0_dp) return
      end if
      keep = .not. ieee_is_nan(x + y)
      a = pack(x, keep)
      b = pack(y, keep)
      r%obs_paired = size(a)
      r%status = mt_insufficient
      if (size(a) < 3) return
      r%status = mt_nonfinite
      if (any(.not. ieee_is_finite(a)) .or. any(.not. ieee_is_finite(b))) return
      if (maxval(abs(a)) > sqrt(huge(1.0_dp)) / real(size(a), dp)) a = a / maxval(abs(a))
      if (maxval(abs(b)) > sqrt(huge(1.0_dp)) / real(size(b), dp)) b = b / maxval(abs(b))
      a = a - r_mean(a)
      b = b - r_mean(b)
      if (maxval(abs(a)) > 0.0_dp) a = a / maxval(abs(a))
      if (maxval(abs(b)) > 0.0_dp) b = b / maxval(abs(b))
      sx = sqrt(sum(a * a))
      sy = sqrt(sum(b * b))
      r%status = mt_degenerate
      if (sx == 0.0_dp .or. sy == 0.0_dp) return
      r%cor = sum((a / sx) * (b / sy))
      if (abs(r%cor - 1.0_dp) < sqrt(epsilon(1.0_dp))) r%cor = 1.0_dp
      if (abs(r%cor + 1.0_dp) < sqrt(epsilon(1.0_dp))) r%cor = -1.0_dp
      r%df = real(size(a) - 2, dp)
      if (abs(r%cor) == 1.0_dp) then
         r%statistic = sign(ieee_value(0.0_dp, ieee_positive_inf), r%cor)
      else
         r%statistic = sqrt(r%df) * r%cor / sqrt(1.0_dp - r%cor ** 2)
      end if
      r%pvalue = student_tail(r%statistic, r%df, alternative)
      if (size(a) > 3 .and. .not. ieee_is_nan(conf_level)) then
         z = atanh(r%cor)
         se = 1.0_dp / sqrt(r%df - 1.0_dp)
         if (alternative == 'two.sided') then
            q = r_qnorm((1.0_dp - conf_level) / 2.0_dp, lower_tail = .false.)
            r%conf_low = tanh(z - se * q)
            r%conf_high = tanh(z + se * q)
         else if (alternative == 'less') then
            r%conf_low = -1.0_dp
            r%conf_high = tanh(z + se * r_qnorm(conf_level))
         else
            r%conf_low = tanh(z - se * r_qnorm(conf_level))
            r%conf_high = 1.0_dp
         end if
      end if
      r%status = mt_ok
   end function correlation_kernel

   pure function group_kernel(x, g, mode) result(r)
      real(dp), intent(in) :: x(:) !! Observations, with IEEE NaNs omitted independently in each group.
      integer, intent(in) :: g(:) !! Group labels of length size(x); nonpositive labels omit observations.
      integer, intent(in) :: mode !! 1 ANOVA, 2 Welch, 3 Kruskal-Wallis, 4 Waerden, 5 Bartlett, 6 Fligner, 7 Levene, 8 Brown-Forsythe.
      type(test_result) :: r
      real(dp), allocatable :: a(:), part(:), counts(:), means(:), vars(:), ranks(:), w(:), med(:)
      integer, allocatable :: groups(:), labels(:)
      logical, allocatable :: keep(:)
      real(dp) :: n, k, center, between, within, tmp, variance, ties, factor
      integer :: i, j, ng
      r = empty_result()
      if (size(x) /= size(g)) return
      keep = .not. ieee_is_nan(x) .and. g > 0
      if (mode == 7 .or. mode == 8) keep = keep .and. ieee_is_finite(x)
      a = pack(x, keep)
      groups = pack(g, keep)
      allocate(labels(size(a)), counts(size(a)), means(size(a)), vars(size(a)))
      ng = 0
      do i = 1, size(a)
         if (any(labels(1:ng) == groups(i))) cycle
         ng = ng + 1
         labels(ng) = groups(i)
      end do
      do i = 1, ng
         part = pack(a, groups == labels(i))
         counts(i) = real(size(part), dp)
         means(i) = r_mean(part)
         vars(i) = r_variance(part)
      end do
      if (mode == 2 .or. mode == 5) then
         ! Rebuild directly from surviving labels; avoid changing sample alignment.
         do i = ng, 1, -1
            if (counts(i) >= 2.0_dp) cycle
            means(i:ng - 1) = means(i + 1:ng)
            vars(i:ng - 1) = vars(i + 1:ng)
            counts(i:ng - 1) = counts(i + 1:ng)
            labels(i:ng - 1) = labels(i + 1:ng)
            ng = ng - 1
         end do
         keep = [(any(labels(1:ng) == groups(j)), j = 1, size(a))]
         a = pack(a, keep)
         groups = pack(groups, keep)
      end if
      r%obs_groups = ng
      r%obs_tot = nint(sum(counts(1:ng)))
      n = real(r%obs_tot, dp)
      k = real(ng, dp)
      r%status = mt_insufficient
      if (ng < 2) return
      if (mode /= 3 .and. mode /= 4) then
         r%status = mt_nonfinite
         if (any(.not. ieee_is_finite(a))) return
      end if
      r%status = mt_degenerate
      select case (mode)
      case (2)
         if (any(vars(1:ng) <= 0.0_dp)) return
         w = counts(1:ng) / vars(1:ng)
         center = sum(w * means(1:ng)) / sum(w)
         between = sum(w * (means(1:ng) - center) ** 2)
         tmp = sum((1.0_dp - w / sum(w)) ** 2 / (counts(1:ng) - 1.0_dp)) / (k * k - 1.0_dp)
         r%df_num = k - 1.0_dp
         r%df_denom = 1.0_dp / (3.0_dp * tmp)
         r%statistic = between / ((k - 1.0_dp) * (1.0_dp + 2.0_dp * (k - 2.0_dp) * tmp))
         r%pvalue = r_pf(r%statistic, r%df_num, r%df_denom, lower_tail = .false.)
      case (5)
         r%var_pooled = sum((counts(1:ng) - 1.0_dp) * vars(1:ng)) / (n - k)
         if (r%var_pooled == 0.0_dp) return
         r%df = k - 1.0_dp
         if (any(vars(1:ng) == 0.0_dp)) then
            r%statistic = ieee_value(0.0_dp, ieee_positive_inf)
            r%pvalue = 0.0_dp
         else
            tmp = 1.0_dp + (sum(1.0_dp / (counts(1:ng) - 1.0_dp)) - 1.0_dp / (n - k)) / (3.0_dp * r%df)
            r%statistic = ((n - k) * log(r%var_pooled) - &
               sum((counts(1:ng) - 1.0_dp) * log(vars(1:ng)))) / tmp
            r%pvalue = r_pchisq(r%statistic, r%df, lower_tail = .false.)
         end if
      case (3, 4, 6)
         if (mode == 6) then
            if (maxval(counts(1:ng)) < 2.0_dp) return
            do i = 1, ng
               part = pack(a, groups == labels(i))
               med = row_medians(reshape(part, [1, size(part)]))
               where (groups == labels(i)) a = abs(a - med(1))
            end do
            if (all(a == 0.0_dp)) return
         end if
         call r_average_ranks(a, ranks)
         if (mode == 3) then
            ties = tie_correction(ranks)
            factor = 1.0_dp - ties / (n ** 3 - n)
            if (factor <= 0.0_dp) return
            tmp = 0.0_dp
            do i = 1, ng
               tmp = tmp + sum(ranks, mask = groups == labels(i)) ** 2 / counts(i)
            end do
            r%statistic = (12.0_dp * tmp / (n * (n + 1.0_dp)) - 3.0_dp * (n + 1.0_dp)) / factor
         else
            if (mode == 4) then
               ranks = r_qnorm(ranks / (n + 1.0_dp))
            else
               ranks = r_qnorm((1.0_dp + ranks / (n + 1.0_dp)) / 2.0_dp)
               ranks = ranks - r_mean(ranks)
            end if
            variance = sum(ranks ** 2) / (n - 1.0_dp)
            if (variance == 0.0_dp) return
            tmp = 0.0_dp
            do i = 1, ng
               tmp = tmp + sum(ranks, mask = groups == labels(i)) ** 2 / counts(i)
            end do
            r%statistic = tmp / variance
         end if
         r%df = k - 1.0_dp
         r%pvalue = r_pchisq(r%statistic, r%df, lower_tail = .false.)
      case (1, 7, 8)
         if (mode == 7 .or. mode == 8) then
            if (maxval(counts(1:ng)) < 3.0_dp) return
            do i = 1, ng
               part = pack(a, groups == labels(i))
               if (mode == 7) then
                  center = means(i)
               else
                  med = row_medians(reshape(part, [1, size(part)]))
                  center = med(1)
               end if
               where (groups == labels(i)) a = abs(a - center)
               part = pack(a, groups == labels(i))
               means(i) = r_mean(part)
               vars(i) = r_variance(part)
            end do
         end if
         center = r_mean(a)
         between = sum(counts(1:ng) * (means(1:ng) - center) ** 2)
         within = 0.0_dp
         do i = 1, ng
            if (counts(i) > 1.0_dp) within = within + (counts(i) - 1.0_dp) * vars(i)
         end do
         r%sumsq_between = between
         r%sumsq_within = within
         r%meansq_between = between / (k - 1.0_dp)
         if (n > k) r%meansq_within = within / (n - k)
         if (n == k) return
         if (within == 0.0_dp .and. (between == 0.0_dp .or. mode /= 1)) return
         r%df_num = k - 1.0_dp
         r%df_denom = n - k
         if (within == 0.0_dp) then
            r%statistic = ieee_value(0.0_dp, ieee_positive_inf)
         else
            r%statistic = r%meansq_between / r%meansq_within
         end if
         r%pvalue = r_pf(r%statistic, r%df_num, r%df_denom, lower_tail = .false.)
      end select
      r%status = mt_ok
      if (ieee_is_nan(r%statistic) .or. ieee_is_nan(r%pvalue)) r%status = mt_nonfinite
   end function group_kernel

   pure real(dp) function tie_correction(ranks) result(ties)
      real(dp), intent(in) :: ranks(:) !! Average ranks, with exact equality marking each tie block.
      integer :: i
      integer, allocatable :: counts(:)
      real(dp) :: n
      allocate(counts(2 * size(ranks)))
      counts = 0
      do i = 1, size(ranks)
         counts(nint(2.0_dp * ranks(i))) = counts(nint(2.0_dp * ranks(i))) + 1
      end do
      ties = 0.0_dp
      do i = 1, size(counts)
         if (counts(i) < 2) cycle
         n = real(counts(i), dp)
         ties = ties + n ** 3 - n
      end do
   end function tie_correction

   pure function wilcoxon_kernel(x, y, mode, null, alternative, exact, correct) result(r)
      real(dp), intent(in) :: x(:) !! First sample, with NaNs omitted.
      real(dp), intent(in) :: y(:) !! Second sample or empty array; paired mode requires length size(x).
      integer, intent(in) :: mode !! 1 signed-rank, 2 two-sample rank sum, or 3 paired signed-rank.
      real(dp), intent(in) :: null !! Finite null location shift applied to x or paired differences.
      character(*), intent(in) :: alternative !! two.sided, less, or greater.
      logical, intent(in), optional :: exact !! Request exact inference; absent selects sample sizes below 50.
      logical, intent(in), optional :: correct !! Continuity correction for normal approximation, default true.
      type(test_result) :: r
      type(rank_distribution) :: distribution
      real(dp), allocatable :: a(:), b(:), ranks(:), scores(:)
      real(dp) :: nx, ny, n, ties, variance, z, correction, stat, lower, upper
      logical :: has_zeroes, use_exact, use_correct
      integer :: i, status
      r = empty_result()
      r%null = null
      r%alternative = alternative
      if (.not. valid_alternative(alternative) .or. .not. ieee_is_finite(null)) return
      a = pack(x, .not. ieee_is_nan(x))
      b = pack(y, .not. ieee_is_nan(y))
      r%obs_x = size(a)
      r%obs_y = size(b)
      r%obs_tot = size(a) + size(b)
      has_zeroes = .false.
      if (mode == 2) then
         a = a - null
         scores = [a, b]
      else
         if (mode == 3) then
            if (size(x) /= size(y)) return
            a = x - y
            a = pack(a, .not. ieee_is_nan(a))
         end if
         a = a - null
         has_zeroes = any(a == 0.0_dp)
         a = pack(a, a /= 0.0_dp)
         scores = abs(a)
         r%obs = size(a)
         r%obs_paired = size(a)
      end if
      nx = real(size(a), dp)
      ny = real(size(b), dp)
      n = nx + ny
      r%status = mt_insufficient
      call r_average_ranks(scores, ranks)
      ties = tie_correction(ranks)
      use_exact = size(a) < 50
      if (mode == 2) use_exact = use_exact .and. size(b) < 50
      if (present(exact)) use_exact = exact
      use_exact = use_exact .and. ties == 0.0_dp .and. .not. has_zeroes
      use_correct = .true.
      if (present(correct)) use_correct = correct
      use_correct = use_correct .and. .not. use_exact
      r%exact = use_exact
      r%corrected = use_correct
      if (size(a) < 1) return
      if (mode == 2 .and. size(b) < 1) return
      if (mode == 2) then
         stat = sum(ranks(1:size(a)))
         r%statistic = stat - nx * (nx + 1.0_dp) / 2.0_dp
         z = r%statistic - nx * ny / 2.0_dp
         variance = nx * ny / 12.0_dp * (n + 1.0_dp - ties / (n * (n - 1.0_dp)))
      else
         stat = sum(ranks, mask = a > 0.0_dp)
         r%statistic = stat
         z = stat - nx * (nx + 1.0_dp) / 4.0_dp
         variance = nx * (nx + 1.0_dp) * (2.0_dp * nx + 1.0_dp) / 24.0_dp - ties / 48.0_dp
      end if
      if (use_exact) then
         scores = [(real(i, dp), i = 1, size(ranks))]
         call permutation_distribution(scores, size(a), distribution, status, paired = mode /= 2, factor = 1.0_dp)
         if (status /= rank_ok) then
            r%status = mt_limit
            return
         end if
         lower = pperm(stat, distribution, 'less')
         upper = pperm(stat, distribution, 'greater')
      else
         correction = 0.0_dp
         if (use_correct) then
            if (alternative == 'two.sided') then
               if (z /= 0.0_dp) correction = sign(0.5_dp, z)
            else if (alternative == 'greater') then
               correction = 0.5_dp
            else
               correction = -0.5_dp
            end if
         end if
         if (variance == 0.0_dp) then
            if (z == correction) then
               r%status = mt_degenerate
               return
            end if
            z = sign(ieee_value(0.0_dp, ieee_positive_inf), z - correction)
         else
            z = (z - correction) / sqrt(variance)
         end if
         lower = r_pnorm(z)
         upper = r_pnorm(z, lower_tail = .false.)
      end if
      select case (alternative)
      case ('two.sided')
         r%pvalue = min(1.0_dp, 2.0_dp * min(lower, upper))
      case ('less')
         r%pvalue = lower
      case ('greater')
         r%pvalue = upper
      end select
      r%status = mt_ok
   end function wilcoxon_kernel

   pure function ad_kernel(x) result(r)
      real(dp), intent(in) :: x(:) !! Sample with NaNs omitted, at least eight nonmissing values required.
      type(test_result) :: r
      type(normality_result) :: test
      r = empty_result()
      test = ad_test(x)
      r%obs = test%n
      r%statistic = test%statistic
      r%pvalue = test%p_value
      r%status = mt_ok
      select case (test%status)
      case (1)
         r%status = mt_insufficient
      case (2)
         r%status = mt_nonfinite
      case (3)
         r%status = mt_degenerate
      end select
   end function ad_kernel

   function jarque_kernel(x) result(r)
      real(dp), intent(in) :: x(:) !! Sample with NaNs omitted, at least two finite values required.
      type(test_result) :: r
      type(moments_test_result) :: test
      real(dp), allocatable :: a(:)
      r = empty_result()
      a = pack(x, .not. ieee_is_nan(x))
      r%obs = size(a)
      if (size(a) > 0 .and. all(ieee_is_finite(a))) then
         if (maxval(abs(a)) > sqrt(huge(1.0_dp)) / real(size(a), dp)) a = a / maxval(abs(a))
         a = a - r_mean(a)
         if (maxval(abs(a)) > 0.0_dp) a = a / maxval(abs(a))
      end if
      test = jarque_test(a)
      r%statistic = test%statistic
      r%pvalue = test%p_value
      r%skewness = skewness(a)
      r%kurtosis = kurtosis(a)
      r%status = mt_degenerate
      if (test%status == 2) r%status = mt_insufficient
      if (test%status == 4) r%status = mt_nonfinite
      if (test%status == 0) then
         r%df = 2.0_dp
         r%status = mt_ok
      end if
   end function jarque_kernel

   pure function cosinor_kernel(x, t, period) result(r)
      real(dp), intent(in) :: x(:) !! Measurements; NaNs and infinities are removed with matching times.
      real(dp), intent(in) :: t(:) !! Measurement times in the same units as period; NaN times are omitted.
      real(dp), intent(in) :: period !! Positive finite period, conventionally 24 for hourly times.
      type(test_result) :: r
      real(dp), allocatable :: a(:), times(:), design(:, :), residual(:)
      logical, allocatable :: keep(:)
      real(dp) :: coefficients(3), sstot, ssres, pi
      integer :: rank, status
      r = empty_result()
      r%period = period
      if (size(x) /= size(t)) return
      if (.not. ieee_is_finite(period) .or. period <= 0.0_dp) return
      keep = ieee_is_finite(x) .and. .not. ieee_is_nan(t)
      a = pack(x, keep)
      times = pack(t, keep)
      r%obs = size(a)
      r%status = mt_insufficient
      if (size(a) == 0) return
      r%status = mt_nonfinite
      if (any(.not. ieee_is_finite(times))) return
      pi = acos(-1.0_dp)
      allocate(design(size(a), 3))
      design(:, 1) = 1.0_dp
      design(:, 2) = sin(2.0_dp * pi * modulo(times, period) / period)
      design(:, 3) = cos(2.0_dp * pi * modulo(times, period) / period)
      call fit_cosinor(design, a, coefficients, rank, status)
      if (status /= 0) return
      where (abs(coefficients) < epsilon(1.0_dp)) coefficients = 0.0_dp
      r%mesor = coefficients(1)
      r%amplitude = sqrt(coefficients(2) ** 2 + coefficients(3) ** 2)
      r%acrophase = modulo(atan2(coefficients(2), coefficients(3)) / (2.0_dp * pi) * period + period, period)
      residual = a - matmul(design, coefficients)
      sstot = sum((a - r_mean(a)) ** 2)
      ssres = sum(residual ** 2)
      r%status = mt_nonfinite
      if (.not. ieee_is_finite(sstot) .or. .not. ieee_is_finite(ssres)) return
      if (sstot < epsilon(1.0_dp)) sstot = 0.0_dp
      if (ssres < epsilon(1.0_dp)) ssres = 0.0_dp
      if (abs(ssres - sstot) < sqrt(epsilon(1.0_dp))) ssres = sstot
      r%status = mt_degenerate
      if (sstot == 0.0_dp) return
      r%rsquared = 1.0_dp - ssres / sstot
      if (size(a) <= 3 .or. rank <= 1) return
      r%df_num = real(rank - 1, dp)
      r%df_denom = real(size(a) - rank, dp)
      if (ssres == 0.0_dp) then
         r%statistic = ieee_value(0.0_dp, ieee_positive_inf)
      else
         r%statistic = ((sstot - ssres) / r%df_num) / (ssres / r%df_denom)
      end if
      r%pvalue = r_pf(r%statistic, r%df_num, r%df_denom, lower_tail = .false.)
      r%status = mt_ok
      if (ieee_is_nan(r%statistic) .or. ieee_is_nan(r%pvalue)) r%status = mt_nonfinite
   end function cosinor_kernel
   pure subroutine fit_cosinor(design, observations, coefficients, rank, status)
      real(dp), intent(in) :: design(:, :) !! Cosinor design matrix with three columns.
      real(dp), intent(in) :: observations(:) !! Finite responses, one per design row.
      real(dp), intent(out) :: coefficients(3) !! Minimum-norm intercept, sine, and cosine coefficients.
      integer, intent(out) :: rank !! Numerical design rank at relative tolerance 1e-7.
      integer, intent(out) :: status !! LAPACK status, zero indicates successful least-squares solution.
      real(dp), allocatable :: matrix(:, :), rhs(:, :), work(:)
      real(dp) :: query(1)
      integer :: pivots(3), m, ld, lwork
      m = size(observations)
      ld = max(3, m)
      allocate(matrix(max(1, m), 3), rhs(ld, 1))
      matrix = design
      rhs = 0.0_dp
      rhs(1:m, 1) = observations
      pivots = 0
      call la_dgelsy(m, 3, 1, matrix, max(1, m), rhs, ld, pivots, 1.0e-7_dp, rank, query, -1, status)
      if (status /= 0) return
      lwork = max(1, ceiling(query(1)))
      allocate(work(lwork))
      matrix = design
      rhs = 0.0_dp
      rhs(1:m, 1) = observations
      pivots = 0
      call la_dgelsy(m, 3, 1, matrix, max(1, m), rhs, ld, pivots, 1.0e-7_dp, rank, work, lwork, status)
      coefficients = rhs(1:3, 1)
   end subroutine fit_cosinor
end module matrixtests_kernels
