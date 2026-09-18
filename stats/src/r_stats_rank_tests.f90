! SPDX-License-Identifier: MIT
! SPDX-FileComment: Exact and rank-based tests corresponding to selected R stats functions.
module r_stats_rank_tests
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use r_distributions, only: r_pchisq, r_pnorm
   use r_kinds, only: dp, r_pi
   use r_optional, only: optval
   use r_sorting, only: r_average_ranks, r_sort
   use r_special, only: r_log_choose
   use r_stats_types, only: fisher_test_result_t, kruskal_test_result_t, &
                            ks_test_result_t, wilcox_test_result_t
   implicit none
   private

   public :: fisher_test, kruskal_test, ks_test, ks_test_normal, ks_test_two_sample
   public :: wilcox_test, wilcox_test_two_sample

   interface fisher_test
      module procedure fisher_test_integer, fisher_test_real
   end interface fisher_test

   interface ks_test
      module procedure ks_test_normal, ks_test_two_sample
   end interface ks_test

   interface wilcox_test
      module procedure wilcox_test_two_sample
   end interface wilcox_test

contains

   pure function fisher_test_integer(table) result(test_result)
      !! Performs a two-sided Fisher exact test on an integer two-by-two table.
      integer, intent(in) :: table(:, :) !! Nonnegative two-by-two contingency table.
      type(fisher_test_result_t) :: test_result

      if (size(table, 1) /= 2 .or. size(table, 2) /= 2 .or. any(table < 0)) return
      test_result = fisher_test_counts(table(1, 1), table(1, 2), table(2, 1), table(2, 2))
   end function fisher_test_integer

   pure function fisher_test_real(table) result(test_result)
      !! Performs a two-sided Fisher exact test when real inputs are integer-valued counts.
      real(dp), intent(in) :: table(:, :) !! Nonnegative integer-valued two-by-two table.
      type(fisher_test_result_t) :: test_result
      integer :: counts(2, 2)

      if (size(table, 1) /= 2 .or. size(table, 2) /= 2) return
      if (any(.not. ieee_is_finite(table)) .or. any(table < 0.0_dp)) return
      counts = nint(table)
      if (any(abs(table - real(counts, dp)) > sqrt(epsilon(1.0_dp)))) return
      test_result = fisher_test_integer(counts)
   end function fisher_test_real

   pure function fisher_test_counts(a, b, c, d) result(test_result)
      !! Evaluates the two-sided hypergeometric tail and conditional odds-ratio estimate.
      integer, intent(in) :: a !! Upper-left table count.
      integer, intent(in) :: b !! Upper-right table count.
      integer, intent(in) :: c !! Lower-left table count.
      integer, intent(in) :: d !! Lower-right table count.
      type(fisher_test_result_t) :: test_result
      integer :: candidate, column_one, high, low, row_one, row_two, total
      real(dp) :: candidate_log_probability, observed_log_probability

      row_one = a + b
      row_two = c + d
      column_one = a + c
      total = row_one + row_two
      if (total <= 0) return
      low = max(0, column_one - row_two)
      high = min(row_one, column_one)
      observed_log_probability = hypergeometric_log_probability(a, row_one, row_two, column_one, total)
      test_result%p_value = 0.0_dp
      do candidate = low, high
         candidate_log_probability = hypergeometric_log_probability(candidate, row_one, row_two, &
                                                                      column_one, total)
         if (candidate_log_probability <= observed_log_probability + &
             1000.0_dp * epsilon(1.0_dp)) then
            test_result%p_value = test_result%p_value + exp(candidate_log_probability)
         end if
      end do
      test_result%p_value = min(test_result%p_value, 1.0_dp)
      test_result%estimate = conditional_odds_ratio(a, row_one, row_two, column_one)
   end function fisher_test_counts

   pure function hypergeometric_log_probability(cell, row_one, row_two, column_one, total) &
      result(log_probability)
      !! Returns a central hypergeometric log probability for fixed table margins.
      integer, intent(in) :: cell !! Candidate upper-left count.
      integer, intent(in) :: row_one !! First row total.
      integer, intent(in) :: row_two !! Second row total.
      integer, intent(in) :: column_one !! First column total.
      integer, intent(in) :: total !! Total table count.
      real(dp) :: log_probability

      log_probability = r_log_choose(row_one, cell) + &
                        r_log_choose(row_two, column_one - cell) - &
                        r_log_choose(total, column_one)
   end function hypergeometric_log_probability

   pure function noncentral_hypergeometric_mean(log_odds, row_one, row_two, column_one) &
      result(mean_cell)
      !! Evaluates the noncentral hypergeometric mean by normalized log-scale summation.
      real(dp), intent(in) :: log_odds !! Log odds ratio.
      integer, intent(in) :: row_one !! First row total.
      integer, intent(in) :: row_two !! Second row total.
      integer, intent(in) :: column_one !! First column total.
      real(dp) :: mean_cell
      integer :: candidate, high, low
      real(dp) :: log_weight, maximum_log_weight, weight, weight_sum

      low = max(0, column_one - row_two)
      high = min(row_one, column_one)
      maximum_log_weight = -huge(1.0_dp)
      do candidate = low, high
         log_weight = r_log_choose(row_one, candidate) + &
                      r_log_choose(row_two, column_one - candidate) + &
                      real(candidate, dp) * log_odds
         maximum_log_weight = max(maximum_log_weight, log_weight)
      end do
      mean_cell = 0.0_dp
      weight_sum = 0.0_dp
      do candidate = low, high
         log_weight = r_log_choose(row_one, candidate) + &
                      r_log_choose(row_two, column_one - candidate) + &
                      real(candidate, dp) * log_odds
         weight = exp(log_weight - maximum_log_weight)
         mean_cell = mean_cell + real(candidate, dp) * weight
         weight_sum = weight_sum + weight
      end do
      mean_cell = mean_cell / weight_sum
   end function noncentral_hypergeometric_mean

   pure function conditional_odds_ratio(observed, row_one, row_two, column_one) result(odds_ratio)
      !! Finds the conditional maximum-likelihood odds ratio by monotone bisection.
      integer, intent(in) :: observed !! Observed upper-left count.
      integer, intent(in) :: row_one !! First row total.
      integer, intent(in) :: row_two !! Second row total.
      integer, intent(in) :: column_one !! First column total.
      real(dp) :: odds_ratio
      integer :: high_count, iteration, low_count
      real(dp) :: high_log_odds, low_log_odds, midpoint

      low_count = max(0, column_one - row_two)
      high_count = min(row_one, column_one)
      if (observed <= low_count) then
         odds_ratio = 0.0_dp
         return
      end if
      if (observed >= high_count) then
         odds_ratio = huge(1.0_dp)
         return
      end if
      low_log_odds = -60.0_dp
      high_log_odds = 60.0_dp
      do iteration = 1, 120
         midpoint = 0.5_dp * (low_log_odds + high_log_odds)
         if (noncentral_hypergeometric_mean(midpoint, row_one, row_two, column_one) < &
             real(observed, dp)) then
            low_log_odds = midpoint
         else
            high_log_odds = midpoint
         end if
      end do
      odds_ratio = exp(0.5_dp * (low_log_odds + high_log_odds))
   end function conditional_odds_ratio

   pure function wilcox_test_two_sample(x, y, paired, correct) result(test_result)
      !! Performs an asymptotic Wilcoxon rank-sum or paired signed-rank test.
      real(dp), intent(in) :: x(:) !! First sample or first paired measurements.
      real(dp), intent(in) :: y(:) !! Second sample or second paired measurements.
      logical, intent(in), optional :: paired !! Use paired signed ranks when true.
      logical, intent(in), optional :: correct !! Apply continuity correction; defaults to true.
      type(wilcox_test_result_t) :: test_result
      real(dp), allocatable :: absolute_differences(:), combined(:), differences(:), ranks(:)
      real(dp) :: center, correction, deviation, standard_deviation, tie_sum
      integer :: first_count, i, nonzero_count, second_count

      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
      if (optval(paired, .false.)) then
         if (size(x) /= size(y) .or. size(x) == 0) return
         differences = x - y
         nonzero_count = count(differences /= 0.0_dp)
         if (nonzero_count == 0) return
         absolute_differences = pack(abs(differences), differences /= 0.0_dp)
         call r_average_ranks(absolute_differences, ranks)
         test_result%statistic = 0.0_dp
         first_count = 0
         do i = 1, size(differences)
            if (differences(i) == 0.0_dp) cycle
            first_count = first_count + 1
            if (differences(i) > 0.0_dp) test_result%statistic = test_result%statistic + ranks(first_count)
         end do
         center = real(nonzero_count, dp) * real(nonzero_count + 1, dp) / 4.0_dp
         tie_sum = rank_tie_sum(absolute_differences)
         standard_deviation = sqrt((real(nonzero_count, dp) * real(nonzero_count + 1, dp) * &
                                   real(2 * nonzero_count + 1, dp) - 0.5_dp * tie_sum) / 24.0_dp)
         test_result%method = 2
      else
         first_count = size(x)
         second_count = size(y)
         if (first_count == 0 .or. second_count == 0) return
         combined = [x, y]
         call r_average_ranks(combined, ranks)
         test_result%statistic = sum(ranks(:first_count)) - &
                                 real(first_count * (first_count + 1), dp) / 2.0_dp
         center = real(first_count, dp) * real(second_count, dp) / 2.0_dp
         tie_sum = rank_tie_sum(combined)
         standard_deviation = sqrt(real(first_count, dp) * real(second_count, dp) / 12.0_dp * &
            (real(first_count + second_count + 1, dp) - tie_sum / &
             (real(first_count + second_count, dp) * &
              real(first_count + second_count - 1, dp))))
         test_result%method = 1
      end if
      if (standard_deviation <= 0.0_dp) return
      deviation = test_result%statistic - center
      correction = 0.0_dp
      if (optval(correct, .true.) .and. deviation /= 0.0_dp) correction = 0.5_dp * sign(1.0_dp, deviation)
      test_result%p_value = 2.0_dp * r_pnorm(abs((deviation - correction) / standard_deviation), &
                                             lower_tail=.false.)
   end function wilcox_test_two_sample

   pure function rank_tie_sum(values) result(tie_sum)
      !! Returns the sum of `t**3-t` over tied groups of ranks.
      real(dp), intent(in) :: values(:) !! Finite values whose exact ties are counted.
      real(dp) :: tie_sum
      real(dp), allocatable :: sorted(:)
      integer :: left, right, tie_count

      call r_sort(values, sorted)
      tie_sum = 0.0_dp
      left = 1
      do while (left <= size(sorted))
         right = left
         do while (right < size(sorted))
            if (sorted(right + 1) /= sorted(left)) exit
            right = right + 1
         end do
         tie_count = right - left + 1
         tie_sum = tie_sum + real(tie_count, dp)**3 - real(tie_count, dp)
         left = right + 1
      end do
   end function rank_tie_sum

   pure function kruskal_test(x, groups) result(test_result)
      !! Performs a tie-corrected Kruskal-Wallis test for integer-coded groups.
      real(dp), intent(in) :: x(:) !! Sample observations.
      integer, intent(in) :: groups(:) !! Group labels conforming with `x`.
      type(kruskal_test_result_t) :: test_result
      integer, allocatable :: group_counts(:), unique_groups(:)
      real(dp), allocatable :: rank_sums(:), ranks(:)
      integer :: group_count, group_index, i, observation_count
      real(dp) :: correction, tie_sum

      observation_count = size(x)
      if (size(groups) /= observation_count .or. observation_count <= 1) return
      if (any(.not. ieee_is_finite(x))) return
      allocate (unique_groups(observation_count))
      group_count = 0
      do i = 1, observation_count
         if (.not. any(unique_groups(:group_count) == groups(i))) then
            group_count = group_count + 1
            unique_groups(group_count) = groups(i)
         end if
      end do
      if (group_count <= 1) return
      allocate (group_counts(group_count), source=0)
      allocate (rank_sums(group_count), source=0.0_dp)
      call r_average_ranks(x, ranks)
      do i = 1, observation_count
         do group_index = 1, group_count
            if (groups(i) /= unique_groups(group_index)) cycle
            group_counts(group_index) = group_counts(group_index) + 1
            rank_sums(group_index) = rank_sums(group_index) + ranks(i)
            exit
         end do
      end do
      test_result%statistic = 12.0_dp / real(observation_count * (observation_count + 1), dp) * &
                              sum(rank_sums**2 / real(group_counts, dp)) - &
                              3.0_dp * real(observation_count + 1, dp)
      tie_sum = rank_tie_sum(x)
      correction = 1.0_dp - tie_sum / &
                   (real(observation_count, dp)**3 - real(observation_count, dp))
      if (correction <= 0.0_dp) return
      test_result%statistic = test_result%statistic / correction
      test_result%parameter = group_count - 1
      test_result%p_value = r_pchisq(test_result%statistic, real(test_result%parameter, dp), &
                                     lower_tail=.false.)
   end function kruskal_test

   pure function ks_test_normal(x, mean, sd) result(test_result)
      !! Performs an asymptotic one-sample Kolmogorov-Smirnov test against a normal distribution.
      real(dp), intent(in) :: x(:) !! Sample observations.
      real(dp), intent(in), optional :: mean !! Reference normal mean; defaults to zero.
      real(dp), intent(in), optional :: sd !! Positive reference normal standard deviation.
      type(ks_test_result_t) :: test_result
      real(dp), allocatable :: sorted(:)
      real(dp) :: cdf_value, location, scale
      integer :: i, sample_size

      sample_size = size(x)
      test_result%n = sample_size
      if (sample_size == 0 .or. any(.not. ieee_is_finite(x))) return
      location = optval(mean, 0.0_dp)
      scale = optval(sd, 1.0_dp)
      if (.not. ieee_is_finite(location) .or. .not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
      call r_sort(x, sorted)
      do i = 1, sample_size
         cdf_value = r_pnorm(sorted(i), location=location, scale=scale)
         test_result%statistic = max(test_result%statistic, &
            real(i, dp) / real(sample_size, dp) - cdf_value, &
            cdf_value - real(i - 1, dp) / real(sample_size, dp))
      end do
      test_result%p_value = kolmogorov_survival(sqrt(real(sample_size, dp)) * test_result%statistic)
      test_result%method = 1
   end function ks_test_normal

   pure function ks_test_two_sample(x, y) result(test_result)
      !! Performs an asymptotic two-sample Kolmogorov-Smirnov test.
      real(dp), intent(in) :: x(:) !! First sample observations.
      real(dp), intent(in) :: y(:) !! Second sample observations.
      type(ks_test_result_t) :: test_result
      real(dp), allocatable :: pooled(:), sorted(:)
      real(dp) :: effective_size
      integer :: first_count, i, second_count

      first_count = size(x)
      second_count = size(y)
      test_result%n = first_count + second_count
      test_result%method = 2
      if (first_count == 0 .or. second_count == 0) return
      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
      pooled = [x, y]
      call r_sort(pooled, sorted)
      do i = 1, size(sorted)
         test_result%statistic = max(test_result%statistic, &
            abs(real(count(x <= sorted(i)), dp) / real(first_count, dp) - &
                real(count(y <= sorted(i)), dp) / real(second_count, dp)))
      end do
      effective_size = sqrt(real(first_count * second_count, dp) / &
                            real(first_count + second_count, dp))
      test_result%p_value = kolmogorov_survival(effective_size * test_result%statistic)
   end function ks_test_two_sample

   pure function kolmogorov_survival(z) result(probability)
      !! Evaluates the limiting two-sided Kolmogorov survival function by its alternating series.
      real(dp), intent(in) :: z !! Nonnegative scaled Kolmogorov statistic.
      real(dp) :: distribution_value, probability, term
      integer :: series_index

      if (z <= 0.0_dp) then
         probability = 1.0_dp
         return
      end if
      if (z <= 0.82_dp) then
         distribution_value = 0.0_dp
         do series_index = 1, 100000
            term = exp(-real(2 * series_index - 1, dp)**2 * r_pi**2 / (8.0_dp * z**2))
            distribution_value = distribution_value + term
            if (term <= epsilon(1.0_dp)) exit
         end do
         distribution_value = sqrt(2.0_dp * r_pi) * distribution_value / z
         probability = max(0.0_dp, min(1.0_dp, 1.0_dp - distribution_value))
         return
      end if
      probability = 0.0_dp
      do series_index = 1, 100000
         term = 2.0_dp * exp(-2.0_dp * real(series_index, dp)**2 * z**2)
         if (mod(series_index, 2) == 1) then
            probability = probability + term
         else
            probability = probability - term
         end if
         if (term <= epsilon(1.0_dp)) exit
      end do
      probability = max(0.0_dp, min(1.0_dp, probability))
   end function kolmogorov_survival

end module r_stats_rank_tests
