! SPDX-License-Identifier: MIT
! SPDX-FileComment: Formula-free analysis-of-variance procedures corresponding to R stats.
module r_stats_anova
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_distributions, only: r_pf
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_stats_types, only: anova_comparison_t, lm_fit_t, one_way_anova_t
   implicit none
   private

   public :: anova_lm, one_way_anova
   public :: oneway_test, oneway_test_result_t

   type :: oneway_test_result_t
      !! Stores Welch or pooled-variance one-way inference with real-valued degrees of freedom.
      real(dp) :: statistic = 0.0_dp !! F statistic.
      real(dp) :: df_numerator = 0.0_dp !! Number of groups minus one.
      real(dp) :: df_denominator = 0.0_dp !! Residual or Welch approximate degrees of freedom.
      real(dp) :: p_value = 1.0_dp !! Upper-tail F probability.
      logical :: var_equal = .false. !! Whether pooled-variance inference was requested.
      integer :: status = 0 !! Zero on success, one for invalid input, two for undefined variance.
   end type oneway_test_result_t

contains

   pure function oneway_test(values, groups, var_equal) result(test)
      !! Tests equality of group means, using Welch's unequal-variance method by default.
      real(dp), intent(in) :: values(:) !! Finite observations with at least two per group.
      integer, intent(in) :: groups(:) !! Arbitrary integer group labels, one per observation.
      logical, intent(in), optional :: var_equal !! Select pooled-variance ANOVA when true.
      type(oneway_test_result_t) :: test !! F statistic, real degrees of freedom, and probability.
      type(one_way_anova_t) :: ordinary
      integer, allocatable :: labels(:), membership(:)
      real(dp), allocatable :: variances(:), weights(:)
      real(dp) :: weighted_center, correction, total_weight
      integer :: n, k, i, j

      test%var_equal = optval(var_equal, .false.)
      n = size(values)
      if (size(groups) /= n .or. n < 4 .or. .not. all(ieee_is_finite(values))) then
         test%status = 1
         return
      end if
      allocate (labels(n), membership(n))
      k = 0
      do i = 1, n
         j = 1
         do while (j <= k)
            if (labels(j) == groups(i)) exit
            j = j + 1
         end do
         if (j > k) then
            k = j
            labels(j) = groups(i)
         end if
         membership(i) = j
      end do
      ordinary = one_way_anova(values, membership)
      if (ordinary%status /= 0) then
         test%status = 1
         return
      end if
      if (any(ordinary%group_counts < 2)) then
         test%status = 1
         return
      end if
      test%df_numerator = real(k - 1, dp)
      if (test%var_equal) then
         test%df_denominator = real(n - k, dp)
         test%statistic = ordinary%f_statistic
         test%p_value = ordinary%p_value
         if (ordinary%ss_within == 0.0_dp .and. ordinary%ss_between == 0.0_dp) then
            test%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
            test%p_value = test%statistic
            test%status = 2
         end if
         return
      end if
      allocate (variances(k), source=0.0_dp)
      do i = 1, n
         j = membership(i)
         variances(j) = variances(j) + (values(i) - ordinary%group_means(j))**2
      end do
      variances = variances/real(ordinary%group_counts - 1, dp)
      if (any(variances <= 0.0_dp)) then
         test%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
         test%p_value = test%statistic
         test%df_denominator = test%statistic
         test%status = 2
         return
      end if
      weights = real(ordinary%group_counts, dp)/variances
      total_weight = sum(weights)
      weighted_center = sum(weights*ordinary%group_means)/total_weight
      correction = sum((1.0_dp - weights/total_weight)**2/ &
                       real(ordinary%group_counts - 1, dp))/(real(k, dp)**2 - 1.0_dp)
      test%df_denominator = 1.0_dp/(3.0_dp*correction)
      test%statistic = sum(weights*(ordinary%group_means - weighted_center)**2)/ &
                       (real(k - 1, dp)*(1.0_dp + 2.0_dp*real(k - 2, dp)*correction))
      test%p_value = r_pf(test%statistic, test%df_numerator, test%df_denominator, lower_tail=.false.)
   end function oneway_test

   pure function one_way_anova(values, groups) result(result)
      !! Computes an ordinary one-way ANOVA for contiguous one-based group labels.
      real(dp), intent(in) :: values(:) !! Finite response observations.
      integer, intent(in) :: groups(:) !! Group labels containing every integer from one to `k`.
      type(one_way_anova_t) :: result
      real(dp) :: grand_mean
      integer :: group, i, k, n

      n = size(values)
      if (size(groups) /= n .or. n < 2 .or. .not. all(ieee_is_finite(values)) .or. &
          any(groups < 1)) then
         result%status = 1
         return
      end if
      k = maxval(groups)
      if (k < 2 .or. n <= k) then
         result%status = 1
         return
      end if
      allocate (result%group_means(k), source=0.0_dp)
      allocate (result%group_counts(k), source=0)
      do i = 1, n
         result%group_means(groups(i)) = result%group_means(groups(i)) + values(i)
         result%group_counts(groups(i)) = result%group_counts(groups(i)) + 1
      end do
      if (any(result%group_counts == 0)) then
         result%status = 1
         return
      end if
      result%group_means = result%group_means/real(result%group_counts, dp)
      grand_mean = sum(values)/real(n, dp)
      do group = 1, k
         result%ss_between = result%ss_between + real(result%group_counts(group), dp)* &
                             (result%group_means(group) - grand_mean)**2
      end do
      do i = 1, n
         result%ss_within = result%ss_within + (values(i) - result%group_means(groups(i)))**2
      end do
      result%df_between = k - 1
      result%df_within = n - k
      result%ms_between = result%ss_between/real(result%df_between, dp)
      result%ms_within = result%ss_within/real(result%df_within, dp)
      if (result%ms_within > 0.0_dp) then
         result%f_statistic = result%ms_between/result%ms_within
         result%p_value = r_pf(result%f_statistic, real(result%df_between, dp), &
                               real(result%df_within, dp), lower_tail=.false.)
      else if (result%ms_between > 0.0_dp) then
         result%f_statistic = huge(1.0_dp)
         result%p_value = 0.0_dp
      end if
   end function one_way_anova

   pure function anova_lm(reduced, full) result(result)
      !! Compares nested linear fits using an extra-sum-of-squares F test.
      type(lm_fit_t), intent(in) :: reduced !! Smaller model fitted to the same observations.
      type(lm_fit_t), intent(in) :: full !! Larger model fitted to the same observations.
      type(anova_comparison_t) :: result
      real(dp) :: full_rss, reduced_rss

      if (size(reduced%y) /= size(full%y) .or. reduced%rank >= full%rank .or. &
          full%df <= 0 .or. any(abs(reduced%y - full%y) > 0.0_dp) .or. &
          any(abs(reduced%weights - full%weights) > 0.0_dp)) then
         result%status = 1
         return
      end if
      reduced_rss = sum(reduced%weights*reduced%resid**2)
      full_rss = sum(full%weights*full%resid**2)
      if (reduced_rss < full_rss) then
         result%status = 1
         return
      end if
      result%df_numerator = full%rank - reduced%rank
      result%df_denominator = full%df
      result%sum_of_squares = reduced_rss - full_rss
      result%mean_square = result%sum_of_squares/real(result%df_numerator, dp)
      if (full_rss > 0.0_dp) then
         result%f_statistic = result%mean_square/(full_rss/real(full%df, dp))
         result%p_value = r_pf(result%f_statistic, real(result%df_numerator, dp), &
                               real(result%df_denominator, dp), lower_tail=.false.)
      else if (result%sum_of_squares > 0.0_dp) then
         result%f_statistic = huge(1.0_dp)
         result%p_value = 0.0_dp
      end if
   end function anova_lm

end module r_stats_anova
