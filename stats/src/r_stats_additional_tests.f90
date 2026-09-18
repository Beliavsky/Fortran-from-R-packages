! SPDX-License-Identifier: MIT
! SPDX-FileComment: Additional pure hypothesis tests corresponding to R stats.
module r_stats_additional_tests
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_base_distributions, only: r_dbinom, r_dpois, r_qbeta, r_qpois
   use r_distributions, only: r_pchisq, r_pf, r_qchisq, r_qf
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_sorting, only: r_average_ranks
   use r_stats_types, only: chisq_test_result_t, exact_count_test_result_t, &
                            variance_test_result_t
   implicit none
   private

   public :: bartlett_test, binom_test, friedman_test, poisson_test, var_test
   public :: mcnemar_test

   interface mcnemar_test
      module procedure mcnemar_test_integer, mcnemar_test_real
   end interface mcnemar_test

contains

   pure function mcnemar_test_integer(table, correct) result(test)
      !! Tests symmetry of an integer square contingency table using McNemar's statistic.
      integer, intent(in) :: table(:, :) !! Nonnegative counts with matching row and column categories.
      logical, intent(in), optional :: correct !! Apply the two-category continuity correction; defaults to true.
      type(chisq_test_result_t) :: test !! Statistic, degrees of freedom, and upper-tail probability.

      test = mcnemar_test_real(real(table, dp), correct)
   end function mcnemar_test_integer

   pure function mcnemar_test_real(table, correct) result(test)
      !! Tests symmetry using squared differences between opposite off-diagonal cells.
      real(dp), intent(in) :: table(:, :) !! Finite nonnegative square table with at least two categories.
      logical, intent(in), optional :: correct !! Apply the two-category continuity correction; defaults to true.
      type(chisq_test_result_t) :: test !! Undefined statistics return NaNs when an opposite-cell pair sums to zero.
      real(dp) :: difference, pair_total
      integer :: categories, i, j

      categories = size(table, 1)
      if (categories < 2 .or. size(table, 2) /= categories) error stop 'mcnemar_test: invalid table shape'
      if (.not. all(ieee_is_finite(table)) .or. any(table < 0.0_dp)) then
         error stop 'mcnemar_test: table entries must be finite and nonnegative'
      end if
      test%parameter = categories*(categories - 1)/2
      test%method = 5
      test%statistic = 0.0_dp
      do j = 2, categories
         do i = 1, j - 1
            pair_total = table(i, j) + table(j, i)
            if (pair_total == 0.0_dp) then
               test%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
               test%p_value = test%statistic
               return
            end if
            difference = abs(table(i, j) - table(j, i))
            if (categories == 2 .and. optval(correct, .true.) .and. difference /= 0.0_dp) then
               difference = difference - 1.0_dp
            end if
            test%statistic = test%statistic + difference**2/pair_total
         end do
      end do
      test%p_value = r_pchisq(test%statistic, real(test%parameter, dp), lower_tail=.false.)
   end function mcnemar_test_real

   pure function binom_test(successes, trials, probability, confidence_level) result(test)
      !! Performs an exact two-sided binomial test with a Clopper-Pearson interval.
      integer, intent(in) :: successes !! Observed number of successes.
      integer, intent(in) :: trials !! Number of Bernoulli trials.
      real(dp), intent(in), optional :: probability !! Null success probability; defaults to one half.
      real(dp), intent(in), optional :: confidence_level !! Confidence level; defaults to 0.95.
      type(exact_count_test_result_t) :: test !! Exact test result and confidence interval.
      real(dp) :: alpha, null_probability, observed_mass, threshold
      integer :: count

      null_probability = optval(probability, 0.5_dp)
      alpha = 1.0_dp - optval(confidence_level, 0.95_dp)
      test%status = 1
      if (trials <= 0 .or. successes < 0 .or. successes > trials .or. &
          null_probability <= 0.0_dp .or. null_probability >= 1.0_dp .or. &
          alpha <= 0.0_dp .or. alpha >= 1.0_dp) return
      test%statistic = real(successes, dp)
      test%parameter = real(trials, dp)
      test%estimate = real(successes, dp)/real(trials, dp)
      test%null_value = null_probability
      observed_mass = r_dbinom(successes, trials, null_probability)
      threshold = observed_mass*(1.0_dp + 1.0e-7_dp)
      test%p_value = 0.0_dp
      do count = 0, trials
         if (r_dbinom(count, trials, null_probability) <= threshold) then
            test%p_value = test%p_value + r_dbinom(count, trials, null_probability)
         end if
      end do
      test%p_value = min(1.0_dp, test%p_value)
      if (successes == 0) then
         test%conf_low = 0.0_dp
      else
         test%conf_low = r_qbeta(0.5_dp*alpha, real(successes, dp), &
                                 real(trials - successes + 1, dp))
      end if
      if (successes == trials) then
         test%conf_high = 1.0_dp
      else
         test%conf_high = r_qbeta(1.0_dp - 0.5_dp*alpha, real(successes + 1, dp), &
                                  real(trials - successes, dp))
      end if
      test%method = 1
      test%status = 0
   end function binom_test

   pure function poisson_test(events, exposure, rate, confidence_level) result(test)
      !! Performs an exact two-sided one-sample Poisson rate test.
      integer, intent(in) :: events !! Observed nonnegative event count.
      real(dp), intent(in) :: exposure !! Positive observation exposure.
      real(dp), intent(in), optional :: rate !! Nonnegative null event rate; defaults to one.
      real(dp), intent(in), optional :: confidence_level !! Confidence level; defaults to 0.95.
      type(exact_count_test_result_t) :: test !! Exact test result and confidence interval.
      real(dp) :: alpha, lambda, null_rate, observed_mass, threshold
      integer :: count, upper_count

      null_rate = optval(rate, 1.0_dp)
      alpha = 1.0_dp - optval(confidence_level, 0.95_dp)
      test%status = 1
      if (events < 0 .or. exposure <= 0.0_dp .or. .not. ieee_is_finite(exposure) .or. &
          null_rate < 0.0_dp .or. .not. ieee_is_finite(null_rate) .or. &
          alpha <= 0.0_dp .or. alpha >= 1.0_dp) return
      lambda = exposure*null_rate
      test%statistic = real(events, dp)
      test%parameter = exposure
      test%estimate = real(events, dp)/exposure
      test%null_value = null_rate
      observed_mass = r_dpois(events, lambda)
      threshold = observed_mass*(1.0_dp + 1.0e-7_dp)
      upper_count = int(r_qpois(1.0_dp - epsilon(1.0_dp), lambda))
      upper_count = max(upper_count, events)
      test%p_value = 0.0_dp
      do count = 0, upper_count
         if (r_dpois(count, lambda) <= threshold) then
            test%p_value = test%p_value + r_dpois(count, lambda)
         end if
      end do
      test%p_value = min(1.0_dp, test%p_value)
      if (events == 0) then
         test%conf_low = 0.0_dp
      else
         test%conf_low = 0.5_dp*r_qchisq(0.5_dp*alpha, real(2*events, dp))/exposure
      end if
      test%conf_high = 0.5_dp*r_qchisq(1.0_dp - 0.5_dp*alpha, &
                                      real(2*(events + 1), dp))/exposure
      test%method = 2
      test%status = 0
   end function poisson_test

   pure function var_test(x, y, ratio, confidence_level) result(test)
      !! Performs a two-sided F test for the ratio of two normal-population variances.
      real(dp), intent(in) :: x(:) !! First finite sample with at least two observations.
      real(dp), intent(in) :: y(:) !! Second finite sample with at least two observations.
      real(dp), intent(in), optional :: ratio !! Null variance ratio; defaults to one.
      real(dp), intent(in), optional :: confidence_level !! Confidence level; defaults to 0.95.
      type(variance_test_result_t) :: test !! Variance-ratio test and confidence interval.
      real(dp) :: alpha, lower_probability, null_ratio, variance_x, variance_y

      null_ratio = optval(ratio, 1.0_dp)
      alpha = 1.0_dp - optval(confidence_level, 0.95_dp)
      test%status = 1
      if (size(x) < 2 .or. size(y) < 2 .or. .not. all(ieee_is_finite(x)) .or. &
          .not. all(ieee_is_finite(y)) .or. null_ratio <= 0.0_dp .or. &
          alpha <= 0.0_dp .or. alpha >= 1.0_dp) return
      variance_x = sum((x - sum(x)/real(size(x), dp))**2)/real(size(x) - 1, dp)
      variance_y = sum((y - sum(y)/real(size(y), dp))**2)/real(size(y) - 1, dp)
      if (variance_x <= 0.0_dp .or. variance_y <= 0.0_dp) return
      test%degrees_freedom1 = real(size(x) - 1, dp)
      test%degrees_freedom2 = real(size(y) - 1, dp)
      test%estimate = variance_x/variance_y
      test%null_value = null_ratio
      test%statistic = test%estimate/null_ratio
      lower_probability = r_pf(test%statistic, test%degrees_freedom1, &
                               test%degrees_freedom2)
      test%p_value = min(1.0_dp, 2.0_dp*min(lower_probability, 1.0_dp - lower_probability))
      test%conf_low = test%estimate/r_qf(1.0_dp - 0.5_dp*alpha, &
                                         test%degrees_freedom1, test%degrees_freedom2)
      test%conf_high = test%estimate/r_qf(0.5_dp*alpha, test%degrees_freedom1, &
                                          test%degrees_freedom2)
      test%status = 0
   end function var_test

   pure function bartlett_test(x, group) result(test)
      !! Performs Bartlett's test of equal variances for explicitly grouped observations.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: group(:) !! Positive integer group label for every observation.
      type(chisq_test_result_t) :: test !! Bartlett statistic, degrees of freedom, and p-value.
      integer, allocatable :: counts(:)
      real(dp), allocatable :: means(:), variances(:)
      real(dp) :: correction, pooled, weighted_log
      integer :: group_count, i, label, nonempty

      if (size(x) /= size(group) .or. size(x) < 4 .or. any(group <= 0) .or. &
          .not. all(ieee_is_finite(x))) return
      group_count = maxval(group)
      allocate (counts(group_count), source=0)
      allocate (means(group_count), source=0.0_dp)
      allocate (variances(group_count), source=0.0_dp)
      do i = 1, size(x)
         label = group(i)
         counts(label) = counts(label) + 1
         means(label) = means(label) + x(i)
      end do
      nonempty = count(counts > 0)
      if (nonempty < 2 .or. any(counts > 0 .and. counts < 2)) return
      do label = 1, group_count
         if (counts(label) > 0) means(label) = means(label)/real(counts(label), dp)
      end do
      do i = 1, size(x)
         label = group(i)
         variances(label) = variances(label) + (x(i) - means(label))**2
      end do
      pooled = 0.0_dp
      weighted_log = 0.0_dp
      correction = 0.0_dp
      do label = 1, group_count
         if (counts(label) == 0) cycle
         variances(label) = variances(label)/real(counts(label) - 1, dp)
         if (variances(label) <= 0.0_dp) return
         pooled = pooled + real(counts(label) - 1, dp)*variances(label)
         weighted_log = weighted_log + real(counts(label) - 1, dp)*log(variances(label))
         correction = correction + 1.0_dp/real(counts(label) - 1, dp)
      end do
      pooled = pooled/real(size(x) - nonempty, dp)
      correction = 1.0_dp + (correction - 1.0_dp/real(size(x) - nonempty, dp))/ &
                   real(3*(nonempty - 1), dp)
      test%statistic = (real(size(x) - nonempty, dp)*log(pooled) - weighted_log)/correction
      test%parameter = nonempty - 1
      test%p_value = r_pchisq(test%statistic, real(test%parameter, dp), lower_tail=.false.)
      test%method = 3
   end function bartlett_test

   pure function friedman_test(x) result(test)
      !! Performs the Friedman rank-sum test for a complete block-by-treatment matrix.
      real(dp), intent(in) :: x(:, :) !! Finite matrix with blocks in rows and treatments in columns.
      type(chisq_test_result_t) :: test !! Friedman statistic, degrees of freedom, and p-value.
      real(dp), allocatable :: ranks(:, :), row_ranks(:)
      real(dp) :: tie_sum, tie_correction
      integer :: blocks, count_equal, i, j, k, treatments
      logical, allocatable :: counted(:)

      blocks = size(x, 1)
      treatments = size(x, 2)
      if (blocks < 2 .or. treatments < 2 .or. .not. all(ieee_is_finite(x))) return
      allocate (ranks(blocks, treatments))
      tie_sum = 0.0_dp
      do i = 1, blocks
         call r_average_ranks(x(i, :), row_ranks)
         ranks(i, :) = row_ranks
         allocate (counted(treatments), source=.false.)
         do j = 1, treatments
            if (counted(j)) cycle
            count_equal = 0
            do k = j, treatments
               if (x(i, k) == x(i, j)) then
                  counted(k) = .true.
                  count_equal = count_equal + 1
               end if
            end do
            tie_sum = tie_sum + real(count_equal**3 - count_equal, dp)
         end do
         deallocate (counted)
      end do
      tie_correction = 1.0_dp - tie_sum/ &
         real(blocks*treatments*(treatments**2 - 1), dp)
      if (tie_correction <= 0.0_dp) return
      test%statistic = (12.0_dp/real(blocks*treatments*(treatments + 1), dp)* &
                        sum(sum(ranks, dim=1)**2) - 3.0_dp*real(blocks*(treatments + 1), dp))/ &
                       tie_correction
      test%parameter = treatments - 1
      test%p_value = r_pchisq(test%statistic, real(test%parameter, dp), lower_tail=.false.)
      test%method = 4
   end function friedman_test

end module r_stats_additional_tests
