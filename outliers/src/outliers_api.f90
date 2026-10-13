! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of outliers 0.15 by Lukasz Komsta.
module outliers_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf
   use r_kinds, only : dp
   use r_descriptive, only : r_mean, r_variance, r_sd
   use r_quantiles, only : r_median, r_quantile_type7
   use r_robust, only : r_mad
   use r_sorting, only : r_sort
   use r_distributions, only : r_pnorm, r_qnorm, r_pt, r_qt, r_pchisq, r_qchisq
   use outliers_distributions
   implicit none
   private
   public :: dp, outlier, outlier_mask, rm_outlier, scores
   public :: grubbs_test, dixon_test, chisq_out_test, cochran_test, cochran_groups
   public :: qtable, qdixon, pdixon, qgrubbs, pgrubbs, qcochran, pcochran
   public :: outlier_test_result, outlier_column

   type :: outlier_test_result
      real(dp) :: statistic = 0.0_dp !! G, Q, chi-square, or Cochran C; U for Grubbs type 20.
      real(dp) :: u = 0.0_dp !! Residual-to-total sum-of-squares ratio for Grubbs tests.
      real(dp) :: p_value = 0.0_dp !! Source-compatible p-value; NaN on invalid input.
      real(dp) :: selected(2) = 0.0_dp !! Selected extreme observations; second used for paired Grubbs tests.
      real(dp) :: df = 0.0_dp !! Mean group size for Cochran, as in the R output.
      integer :: n = 0 !! Number of nonmissing observations, or groups for Cochran.
      integer :: variant = 0 !! Applied Grubbs/Dixon variant; 0 for other tests.
      integer :: group = 0 !! First selected group index for Cochran, with all ties available in groups.
      integer, allocatable :: groups(:) !! All selected Cochran group indices, preserving ties.
      logical :: lower = .false. !! True if the lower endpoint or variance is selected.
      integer :: status = 1 !! 0 success, 1 invalid data/domain; NaN statistics accompany errors.
   end type outlier_test_result

   type :: outlier_column
      real(dp), allocatable :: values(:) !! A column after removal; lengths can differ because of ties.
   end type outlier_column

   interface outlier
      module procedure outlier_vector, outlier_matrix
   end interface
   interface outlier_mask
      module procedure outlier_mask_vector, outlier_mask_matrix
   end interface
   interface rm_outlier
      module procedure rm_outlier_vector, rm_outlier_matrix
   end interface
   interface scores
      module procedure scores_vector, scores_matrix
   end interface
   interface cochran_test
      module procedure cochran_variances, cochran_matrix
   end interface
contains
   pure function invalid_result() result(res)
      type(outlier_test_result) :: res
      res%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      res%u = res%statistic
      res%p_value = res%statistic
      res%selected = res%statistic
      res%df = res%statistic
      allocate(res%groups(0))
   end function invalid_result

   pure logical function choose_lower(x, opposite) result(lower)
      real(dp), intent(in) :: x(:) !! Nonempty, finite sample.
      logical, intent(in), optional :: opposite !! Reverse the more-extreme endpoint selection; default false.
      real(dp) :: center
      center = r_mean(x)
      lower = maxval(x) - center < center - minval(x)
      if (present(opposite)) lower = lower .neqv. opposite
   end function choose_lower

   pure real(dp) function outlier_vector(x, opposite) result(value)
      real(dp), intent(in) :: x(:) !! Sample; NaNs are omitted, while any infinity makes the result NaN.
      logical, intent(in), optional :: opposite !! Select the less-extreme endpoint; default false.
      real(dp), allocatable :: clean(:)
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      clean = pack(x, .not. ieee_is_nan(x))
      if (size(clean) == 0) return
      if (any(.not. ieee_is_finite(clean))) return
      if (choose_lower(clean, opposite)) then
         value = minval(clean)
      else
         value = maxval(clean)
      end if
   end function outlier_vector

   pure function outlier_matrix(x, opposite) result(values)
      real(dp), intent(in) :: x(:, :) !! Matrix with observations in rows, independently processed columns.
      logical, intent(in), optional :: opposite !! Select the less-extreme endpoint in each column; default false.
      real(dp) :: values(size(x, 2))
      integer :: j
      do j = 1, size(x, 2)
         values(j) = outlier_vector(x(:, j), opposite)
      end do
   end function outlier_matrix

   pure function outlier_mask_vector(x, opposite) result(mask)
      real(dp), intent(in) :: x(:) !! Sample; missing positions are false in the returned logical mask.
      logical, intent(in), optional :: opposite !! Select the less-extreme endpoint; default false.
      logical :: mask(size(x))
      real(dp) :: value
      value = outlier_vector(x, opposite)
      mask = .false.
      if (ieee_is_nan(value)) return
      mask = x == value
   end function outlier_mask_vector

   pure function outlier_mask_matrix(x, opposite) result(mask)
      real(dp), intent(in) :: x(:, :) !! Observations by columns; output mask has identical shape.
      logical, intent(in), optional :: opposite !! Reverse endpoint selection in each column; default false.
      logical :: mask(size(x, 1), size(x, 2))
      integer :: j
      do j = 1, size(x, 2)
         mask(:, j) = outlier_mask_vector(x(:, j), opposite)
      end do
   end function outlier_mask_matrix

   pure function rm_outlier_vector(x, fill, median, opposite) result(values)
      real(dp), intent(in) :: x(:) !! Sample; all ties at the selected endpoint are removed or replaced.
      logical, intent(in), optional :: fill !! Replace selected values instead of removing them; default false.
      logical, intent(in), optional :: median !! Fill with remaining-data median instead of mean; default false.
      logical, intent(in), optional :: opposite !! Remove the less-extreme endpoint; default false.
      real(dp), allocatable :: values(:), remaining(:)
      logical :: mask(size(x)), replace, use_median
      real(dp) :: replacement
      replace = .false.
      use_median = .false.
      if (present(fill)) replace = fill
      if (present(median)) use_median = median
      mask = outlier_mask_vector(x, opposite)
      remaining = pack(x, .not. mask)
      if (.not. replace) then
         values = remaining
      else
         values = x
         replacement = r_mean(remaining)
         if (use_median) replacement = r_median(remaining)
         where (mask) values = replacement
      end if
   end function rm_outlier_vector

   pure function rm_outlier_matrix(x, fill, median, opposite) result(columns)
      real(dp), intent(in) :: x(:, :) !! Matrix processed by column; return columns may have unequal lengths.
      logical, intent(in), optional :: fill !! Replace selected values instead of removing them; default false.
      logical, intent(in), optional :: median !! Use the remaining median for replacements; default false.
      logical, intent(in), optional :: opposite !! Remove the less-extreme endpoint; default false.
      type(outlier_column) :: columns(size(x, 2))
      integer :: j
      do j = 1, size(x, 2)
         columns(j)%values = rm_outlier_vector(x(:, j), fill, median, opposite)
      end do
   end function rm_outlier_matrix

   pure elemental real(dp) function ratio(numerator, denominator) result(value)
      real(dp), intent(in) :: numerator !! Dividend; zero divided by zero returns quiet NaN.
      real(dp), intent(in) :: denominator !! Divisor; zero with nonzero dividend returns signed infinity.
      if (ieee_is_nan(numerator) .or. ieee_is_nan(denominator)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (denominator == 0.0_dp) then
         if (numerator == 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            value = sign(ieee_value(0.0_dp, ieee_positive_inf), numerator)
         end if
      else
         value = numerator / denominator
      end if
   end function ratio

   pure function scores_vector(x, type, prob, lim) result(values)
      real(dp), intent(in) :: x(:) !! Sample; NaNs propagate to all scores, as in R's default missing-data behavior.
      character(len = *), intent(in), optional :: type !! Exact score name: z (default), t, chisq, iqr, or mad.
      real(dp), intent(in), optional :: prob !! Absent/NaN: scores; 1: CDF; otherwise 0/1 threshold indicators.
      real(dp), intent(in), optional :: lim !! IQR-score absolute threshold; absent/NaN returns IQR scores.
      real(dp) :: values(size(x)), n, center, spread, q1, q3, threshold, p
      character(len = 5) :: kind
      integer :: i
      logical :: probability
      values = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) return
      n = real(size(x), dp)
      kind = 'z'
      if (present(type)) then
         if (len_trim(type) > len(kind)) return
         kind = type
      end if
      select case (kind)
      case ('z', 't', 'chisq')
         if (size(x) < 2) return
         center = r_mean(x)
         spread = r_sd(x)
         values = ratio(x - center, spread)
         if (kind == 'chisq') values = values ** 2
         if (kind == 't') then
            do i = 1, size(x)
               if (n - 1.0_dp - values(i) ** 2 >= 0.0_dp) then
                  values(i) = ratio(values(i) * sqrt(n - 2.0_dp), sqrt(n - 1.0_dp - values(i) ** 2))
               else
                  values(i) = ieee_value(0.0_dp, ieee_quiet_nan)
               end if
            end do
         end if
      case ('mad')
         values = ratio(x - r_median(x), r_mad(x))
      case ('iqr')
         q1 = r_quantile_type7(x, 0.25_dp)
         q3 = r_quantile_type7(x, 0.75_dp)
         values = 0.0_dp
         where (x < q1) values = ratio(x - q1, q3 - q1)
         where (x > q3) values = ratio(x - q3, q3 - q1)
         if (present(lim)) then
            if (.not. ieee_is_nan(lim)) then
               where (.not. ieee_is_nan(values)) values = merge(1.0_dp, 0.0_dp, abs(values) > lim)
            end if
         end if
         return
      case default
         return
      end select
      probability = .false.
      p = 0.0_dp
      if (present(prob)) then
         p = prob
         probability = .not. ieee_is_nan(p)
      end if
      if (.not. probability) return
      if (p < 0.0_dp .or. p > 1.0_dp .or. .not. ieee_is_finite(p)) then
         values = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (p == 1.0_dp) then
         select case (kind)
         case ('z', 'mad')
            values = r_pnorm(values)
         case ('t')
            values = r_pt(values, n - 2.0_dp)
         case ('chisq')
            values = r_pchisq(values, 1.0_dp)
         end select
      else
         if (p == 0.0_dp .and. kind /= 'chisq') then
            threshold = (n - 1.0_dp) / sqrt(n)
         else
            select case (kind)
            case ('z', 'mad')
               threshold = r_qnorm(p)
            case ('t')
               threshold = r_qt(p, n - 2.0_dp)
            case ('chisq')
               threshold = r_qchisq(p, 1.0_dp)
            end select
         end if
         if (ieee_is_nan(threshold)) then
            values = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            where (.not. ieee_is_nan(values)) values = merge(1.0_dp, 0.0_dp, abs(values) > threshold)
         end if
      end if
   end function scores_vector

   pure function scores_matrix(x, type, prob, lim) result(values)
      real(dp), intent(in) :: x(:, :) !! Matrix with observations in rows; each column is scored separately.
      character(len = *), intent(in), optional :: type !! Exact score name: z (default), t, chisq, iqr, or mad.
      real(dp), intent(in), optional :: prob !! Absent/NaN: scores; 1: CDF; otherwise real-valued 0/1 indicators.
      real(dp), intent(in), optional :: lim !! IQR-score threshold; ignored for other score types.
      real(dp) :: values(size(x, 1), size(x, 2))
      integer :: j
      do j = 1, size(x, 2)
         values(:, j) = scores_vector(x(:, j), type, prob, lim)
      end do
   end function scores_matrix

   pure function chisq_out_test(x, variance, opposite) result(res)
      real(dp), intent(in) :: x(:) !! Sample; NaNs are omitted before computing the test and default variance.
      real(dp), intent(in), optional :: variance !! Positive reference variance; defaults to the cleaned sample variance.
      logical, intent(in), optional :: opposite !! Test the less-extreme endpoint; default false.
      type(outlier_test_result) :: res
      real(dp), allocatable :: clean(:)
      real(dp) :: v, extreme
      res = invalid_result()
      clean = pack(x, .not. ieee_is_nan(x))
      res%n = size(clean)
      if (res%n < 1 .or. any(.not. ieee_is_finite(clean))) return
      v = r_variance(clean)
      if (present(variance)) v = variance
      if (v <= 0.0_dp .or. .not. ieee_is_finite(v)) return
      extreme = outlier_vector(clean, opposite)
      res%selected(1) = extreme
      res%lower = choose_lower(clean, opposite)
      res%statistic = (extreme - r_mean(clean)) ** 2 / v
      res%p_value = r_pchisq(res%statistic, 1.0_dp, lower_tail = .false.)
      res%status = 0
   end function chisq_out_test

   pure function dixon_test(x, type, opposite, two_sided) result(res)
      real(dp), intent(in) :: x(:) !! Sample; NaNs are removed before selecting the sample-size-specific ratio.
      integer, intent(in), optional :: type !! Dixon code 0 (default automatic), 10, 11, 12, 20, 21, or 22.
      logical, intent(in), optional :: opposite !! Test the less-extreme endpoint; default false.
      logical, intent(in), optional :: two_sided !! Apply R's doubled-and-reflected p-value; default true.
      type(outlier_test_result) :: res
      real(dp), allocatable :: clean(:)
      real(dp) :: numerator, denominator
      integer :: n, kind, gap, skip
      logical :: two
      res = invalid_result()
      call r_sort(x, clean, na_rm = .true.)
      n = size(clean)
      res%n = n
      kind = 0
      if (present(type)) kind = type
      kind = dixon_type(n, kind)
      if (kind < 0 .or. any(.not. ieee_is_finite(clean))) return
      if (clean(n) <= clean(1)) return
      res%variant = kind
      res%lower = choose_lower(clean, opposite)
      gap = kind / 10
      skip = mod(kind, 10)
      if (res%lower) then
         numerator = clean(1 + gap) - clean(1)
         denominator = clean(n - skip) - clean(1)
         res%selected(1) = clean(1)
      else
         numerator = clean(n) - clean(n - gap)
         denominator = clean(n) - clean(1 + skip)
         res%selected(1) = clean(n)
      end if
      if (denominator <= 0.0_dp) return
      res%statistic = numerator / denominator
      res%p_value = pdixon(res%statistic, n, kind)
      two = .true.
      if (present(two_sided)) two = two_sided
      if (two) res%p_value = folded_probability(res%p_value)
      if (ieee_is_finite(res%p_value)) res%status = 0
   end function dixon_test

   pure elemental real(dp) function folded_probability(p) result(value)
      real(dp), intent(in) :: p !! One-sided p-value in [0,1], doubled and reflected above 1 as in R.
      value = 2.0_dp * p
      if (value > 1.0_dp) value = 2.0_dp - value
   end function folded_probability

   pure function grubbs_test(x, type, opposite, two_sided) result(res)
      real(dp), intent(in) :: x(:) !! Sample; NaNs are omitted and finite observations sorted.
      integer, intent(in), optional :: type !! Grubbs code 10 (default), 11, or 20.
      logical, intent(in), optional :: opposite !! Reverse endpoint selection for types 10 and 20; ignored for 11.
      logical, intent(in), optional :: two_sided !! Apply R's doubled-and-reflected p-value; default false.
      type(outlier_test_result) :: res
      real(dp), allocatable :: clean(:)
      real(dp) :: variance, rn
      integer :: n, kind
      logical :: two
      res = invalid_result()
      call r_sort(x, clean, na_rm = .true.)
      n = size(clean)
      res%n = n
      kind = 10
      if (present(type)) kind = type
      res%variant = kind
      if (n < 3 .or. any(.not. ieee_is_finite(clean))) return
      if (kind == 20 .and. n < 4) return
      if (kind == 20 .and. n > 30) return
      variance = r_variance(clean)
      if (variance <= 0.0_dp .or. .not. ieee_is_finite(variance)) return
      rn = real(n, dp)
      res%lower = choose_lower(clean, opposite)
      select case (kind)
      case (10)
         res%selected(1) = outlier_vector(clean, opposite)
         res%statistic = abs(res%selected(1) - r_mean(clean)) / sqrt(variance)
         if (res%lower) then
            res%u = r_variance(clean(2:)) / variance * (rn - 2.0_dp) / (rn - 1.0_dp)
         else
            res%u = r_variance(clean(:n - 1)) / variance * (rn - 2.0_dp) / (rn - 1.0_dp)
         end if
         res%p_value = 1.0_dp - pgrubbs(res%statistic, n, kind)
      case (11)
         res%selected = [clean(1), clean(n)]
         res%statistic = (clean(n) - clean(1)) / sqrt(variance)
         res%u = r_variance(clean(2:n - 1)) / variance * (rn - 3.0_dp) / (rn - 1.0_dp)
         res%p_value = 1.0_dp - pgrubbs(res%statistic, n, kind)
      case (20)
         if (res%lower) then
            res%selected = clean(:2)
            res%u = r_variance(clean(3:)) / variance * (rn - 3.0_dp) / (rn - 1.0_dp)
         else
            res%selected = clean(n - 1:)
            res%u = r_variance(clean(:n - 2)) / variance * (rn - 3.0_dp) / (rn - 1.0_dp)
         end if
         res%statistic = res%u
         res%p_value = pgrubbs(res%u, n, kind)
      case default
         return
      end select
      two = .false.
      if (present(two_sided)) two = two_sided
      if (two) res%p_value = folded_probability(res%p_value)
      if (ieee_is_finite(res%p_value)) res%status = 0
   end function grubbs_test

   pure function cochran_variances(variances, counts, inlying) result(res)
      real(dp), intent(in) :: variances(:) !! Nonnegative group variances; at least two groups and positive sum.
      real(dp), intent(in) :: counts(:) !! Group sample sizes, same length as variances and each greater than 1.
      logical, intent(in), optional :: inlying !! Select minimum variance instead of maximum; default false.
      type(outlier_test_result) :: res
      real(dp) :: selected
      integer :: k, j
      logical :: inside
      res = invalid_result()
      k = size(variances)
      res%n = k
      if (k < 2 .or. size(counts) /= k) return
      if (any(.not. ieee_is_finite(variances)) .or. any(.not. ieee_is_finite(counts))) return
      if (any(variances < 0.0_dp) .or. sum(variances) <= 0.0_dp .or. any(counts <= 1.0_dp)) return
      inside = .false.
      if (present(inlying)) inside = inlying
      selected = maxval(variances)
      if (inside) selected = minval(variances)
      res%groups = pack([(j, j = 1, k)], variances == selected)
      res%group = res%groups(1)
      res%selected(1) = selected
      res%lower = inside
      res%df = r_mean(counts)
      res%statistic = selected / sum(variances)
      res%p_value = pcochran(res%statistic, res%df, k)
      if (.not. inside) res%p_value = 1.0_dp - res%p_value
      res%status = 0
   end function cochran_variances

   pure function cochran_matrix(x, inlying) result(res)
      real(dp), intent(in) :: x(:, :) !! Equal-size groups in columns, observations in rows; missing data is invalid.
      logical, intent(in), optional :: inlying !! Test the minimum group variance instead of maximum; default false.
      type(outlier_test_result) :: res
      real(dp) :: variances(size(x, 2)), counts(size(x, 2))
      integer :: j
      do j = 1, size(x, 2)
         variances(j) = r_variance(x(:, j))
      end do
      counts = real(size(x, 1), dp)
      res = cochran_variances(variances, counts, inlying)
   end function cochran_matrix

   pure function cochran_groups(x, groups, inlying) result(res)
      real(dp), intent(in) :: x(:) !! Observations for possibly unequal-size groups; missing data is invalid.
      integer, intent(in) :: groups(:) !! Positive consecutive group labels 1 through k, same length as x.
      logical, intent(in), optional :: inlying !! Select the group with minimum variance; default false.
      type(outlier_test_result) :: res
      real(dp), allocatable :: variances(:), counts(:)
      integer :: k, j
      res = invalid_result()
      if (size(x) /= size(groups) .or. size(x) == 0) return
      if (any(groups < 1) .or. any(.not. ieee_is_finite(x))) return
      k = maxval(groups)
      if (k > size(x) / 2) return
      allocate(variances(k), counts(k))
      do j = 1, k
         counts(j) = real(count(groups == j), dp)
         variances(j) = r_variance(pack(x, groups == j))
      end do
      res = cochran_variances(variances, counts, inlying)
   end function cochran_groups
end module outliers_api
