! SPDX-License-Identifier: MIT
! SPDX-FileComment: Classical hypothesis tests corresponding to selected R stats functions.
module r_stats_tests
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_descriptive, only: r_correlation
   use r_distributions, only: r_pchisq, r_pt, r_qnorm, r_qt
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_sorting, only: r_average_ranks
   use r_stats_types, only: chisq_test_result_t, cor_test_result_t, prop_test_result_t, &
                            t_test_result_t
   implicit none
   private

   public :: chisq_test, cor_test, prop_test
   public :: t_test, t_test_one, t_test_p_value, t_test_p_value_one
   public :: t_test_p_value_two, t_test_two

   interface chisq_test
      module procedure chisq_test_int_vector, chisq_test_real_vector
      module procedure chisq_test_int_matrix, chisq_test_real_matrix
   end interface chisq_test

   interface cor_test
      module procedure cor_test_integer, cor_test_real
   end interface cor_test

   interface prop_test
      module procedure prop_test_int_scalar, prop_test_real_scalar
      module procedure prop_test_int_vector, prop_test_real_vector
   end interface prop_test

   interface t_test
      module procedure t_test_one, t_test_two
   end interface t_test

   interface t_test_p_value
      module procedure t_test_p_value_one, t_test_p_value_two
   end interface t_test_p_value

contains

   pure function chisq_test_real_vector(x, p, rescale_p) result(test_result)
      !! Performs a Pearson chi-square goodness-of-fit test for specified category probabilities.
      real(dp), intent(in) :: x(:) !! Nonnegative observed category counts.
      real(dp), intent(in), optional :: p(:) !! Null category probabilities; defaults to equal probabilities.
      logical, intent(in), optional :: rescale_p !! Normalize supplied probabilities when true.
      type(chisq_test_result_t) :: test_result
      real(dp), allocatable :: expected(:), probabilities(:)
      real(dp) :: probability_sum, total
      integer :: category_count

      category_count = size(x)
      if (category_count <= 1 .or. any(x < 0.0_dp)) return
      total = sum(x)
      if (total <= 0.0_dp) return
      if (present(p)) then
         if (size(p) /= category_count .or. any(p < 0.0_dp)) return
         probability_sum = sum(p)
         if (probability_sum <= 0.0_dp) return
         if (.not. optval(rescale_p, .false.) .and. &
             abs(probability_sum - 1.0_dp) > sqrt(epsilon(1.0_dp))) return
         probabilities = p / probability_sum
      else
         allocate (probabilities(category_count), source=1.0_dp / real(category_count, dp))
      end if
      expected = total * probabilities
      if (any(expected <= 0.0_dp)) return

      test_result%statistic = sum((x - expected)**2 / expected)
      test_result%parameter = category_count - 1
      test_result%p_value = r_pchisq(test_result%statistic, real(test_result%parameter, dp), &
                                     lower_tail=.false.)
      test_result%method = 1
   end function chisq_test_real_vector

   pure function chisq_test_int_vector(x, p, rescale_p) result(test_result)
      !! Performs a goodness-of-fit test for integer category counts.
      integer, intent(in) :: x(:) !! Nonnegative observed category counts.
      real(dp), intent(in), optional :: p(:) !! Null category probabilities; defaults to equal probabilities.
      logical, intent(in), optional :: rescale_p !! Normalize supplied probabilities when true.
      type(chisq_test_result_t) :: test_result

      test_result = chisq_test_real_vector(real(x, dp), p, rescale_p)
   end function chisq_test_int_vector

   pure function chisq_test_real_matrix(x, correct) result(test_result)
      !! Performs a Pearson chi-square independence test on a contingency table.
      real(dp), intent(in) :: x(:, :) !! Nonnegative contingency-table counts.
      logical, intent(in), optional :: correct !! Apply Yates correction to a two-by-two table.
      type(chisq_test_result_t) :: test_result
      real(dp), allocatable :: expected(:, :), column_totals(:), row_totals(:)
      real(dp) :: correction, total
      integer :: column_count, row_count

      row_count = size(x, 1)
      column_count = size(x, 2)
      if (row_count <= 1 .or. column_count <= 1 .or. any(x < 0.0_dp)) return
      total = sum(x)
      if (total <= 0.0_dp) return
      row_totals = sum(x, dim=2)
      column_totals = sum(x, dim=1)
      expected = spread(row_totals, 2, column_count) * spread(column_totals, 1, row_count) / total
      if (any(expected <= 0.0_dp)) return

      correction = 0.0_dp
      if (row_count == 2 .and. column_count == 2 .and. optval(correct, .true.)) then
         correction = min(0.5_dp, minval(abs(x - expected)))
      end if
      test_result%statistic = sum((abs(x - expected) - correction)**2 / expected)
      test_result%parameter = (row_count - 1) * (column_count - 1)
      test_result%p_value = r_pchisq(test_result%statistic, real(test_result%parameter, dp), &
                                     lower_tail=.false.)
      test_result%method = 2
   end function chisq_test_real_matrix

   pure function chisq_test_int_matrix(x, correct) result(test_result)
      !! Performs an independence test for an integer contingency table.
      integer, intent(in) :: x(:, :) !! Nonnegative contingency-table counts.
      logical, intent(in), optional :: correct !! Apply Yates correction to a two-by-two table.
      type(chisq_test_result_t) :: test_result

      test_result = chisq_test_real_matrix(real(x, dp), correct)
   end function chisq_test_int_matrix

   pure function prop_test_real_scalar(x, n, p, correct) result(test_result)
      !! Performs a two-sided one-sample Pearson proportion test with a 95 percent score interval.
      real(dp), intent(in) :: x !! Number of successes.
      real(dp), intent(in) :: n !! Number of trials.
      real(dp), intent(in), optional :: p !! Success probability under the null; defaults to one half.
      logical, intent(in), optional :: correct !! Apply continuity correction; defaults to true.
      type(prop_test_result_t) :: test_result
      real(dp) :: adjusted, correction, null_probability, observed, z, z22n

      null_probability = optval(p, 0.5_dp)
      if (n <= 0.0_dp .or. x < 0.0_dp .or. x > n) return
      if (null_probability <= 0.0_dp .or. null_probability >= 1.0_dp) return
      observed = x / n
      correction = 0.0_dp
      if (optval(correct, .true.)) correction = min(0.5_dp, abs(x - n * null_probability))
      test_result%statistic = (abs(x - n * null_probability) - correction)**2 / &
                              (n * null_probability * (1.0_dp - null_probability))
      test_result%parameter = 1
      test_result%p_value = r_pchisq(test_result%statistic, 1.0_dp, lower_tail=.false.)
      test_result%estimate = observed
      test_result%estimates = [observed]
      test_result%null_value = null_probability
      test_result%null_values = [null_probability]
      test_result%method = 1

      z = r_qnorm(0.975_dp)
      z22n = z**2 / (2.0_dp * n)
      adjusted = observed + correction / n
      if (adjusted >= 1.0_dp) then
         test_result%conf_high = 1.0_dp
      else
         test_result%conf_high = (adjusted + z22n + &
            z * sqrt(adjusted * (1.0_dp - adjusted) / n + z22n / (2.0_dp * n))) / &
            (1.0_dp + 2.0_dp * z22n)
      end if
      adjusted = observed - correction / n
      if (adjusted <= 0.0_dp) then
         test_result%conf_low = 0.0_dp
      else
         test_result%conf_low = (adjusted + z22n - &
            z * sqrt(adjusted * (1.0_dp - adjusted) / n + z22n / (2.0_dp * n))) / &
            (1.0_dp + 2.0_dp * z22n)
      end if
   end function prop_test_real_scalar

   pure function prop_test_int_scalar(x, n, p, correct) result(test_result)
      !! Performs a one-sample proportion test for integer success and trial counts.
      integer, intent(in) :: x !! Number of successes.
      integer, intent(in) :: n !! Number of trials.
      real(dp), intent(in), optional :: p !! Success probability under the null; defaults to one half.
      logical, intent(in), optional :: correct !! Apply continuity correction; defaults to true.
      type(prop_test_result_t) :: test_result

      test_result = prop_test_real_scalar(real(x, dp), real(n, dp), p, correct)
   end function prop_test_int_scalar

   pure function prop_test_real_vector(x, n, p, correct) result(test_result)
      !! Performs a two-sided test of equal or specified proportions across samples.
      real(dp), intent(in) :: x(:) !! Success counts for each sample.
      real(dp), intent(in) :: n(:) !! Trial counts conforming with `x`.
      real(dp), intent(in), optional :: p(:) !! Null probabilities conforming with `x`.
      logical, intent(in), optional :: correct !! Apply continuity correction for at most two samples.
      type(prop_test_result_t) :: test_result
      real(dp), allocatable :: expected_failure(:), expected_success(:), probabilities(:)
      real(dp) :: correction, delta, pooled, width, z
      integer :: sample_count

      sample_count = size(x)
      if (sample_count < 1 .or. size(n) /= sample_count) return
      if (any(n <= 0.0_dp) .or. any(x < 0.0_dp) .or. any(x > n)) return
      if (sample_count == 1) then
         if (present(p)) then
            if (size(p) /= 1) return
            test_result = prop_test_real_scalar(x(1), n(1), p(1), correct)
         else
            test_result = prop_test_real_scalar(x(1), n(1), correct=correct)
         end if
         return
      end if

      test_result%estimates = x / n
      test_result%estimate = test_result%estimates(1)
      test_result%estimate2 = test_result%estimates(2)
      correction = 0.0_dp
      if (present(p)) then
         if (size(p) /= sample_count .or. any(p <= 0.0_dp) .or. any(p >= 1.0_dp)) return
         probabilities = p
         test_result%null_values = p
         test_result%null_value = p(1)
         test_result%parameter = sample_count
         if (sample_count <= 2 .and. optval(correct, .true.)) correction = 0.5_dp
      else
         pooled = sum(x) / sum(n)
         if (pooled <= 0.0_dp .or. pooled >= 1.0_dp) return
         allocate (probabilities(sample_count), source=pooled)
         test_result%null_value = pooled
         test_result%parameter = sample_count - 1
         if (sample_count == 2 .and. optval(correct, .true.)) then
            delta = test_result%estimate - test_result%estimate2
            correction = min(0.5_dp, abs(delta) / sum(1.0_dp / n))
         end if
      end if
      expected_success = n * probabilities
      expected_failure = n * (1.0_dp - probabilities)
      test_result%statistic = sum((abs(x - expected_success) - correction)**2 / expected_success)
      test_result%statistic = test_result%statistic + &
         sum((abs(n - x - expected_failure) - correction)**2 / expected_failure)
      test_result%p_value = r_pchisq(test_result%statistic, real(test_result%parameter, dp), &
                                     lower_tail=.false.)
      test_result%method = 2

      if (sample_count == 2 .and. .not. present(p)) then
         delta = test_result%estimate - test_result%estimate2
         z = r_qnorm(0.975_dp)
         width = z * sqrt(sum(test_result%estimates * (1.0_dp - test_result%estimates) / n)) + &
                 correction * sum(1.0_dp / n)
         test_result%conf_low = max(delta - width, -1.0_dp)
         test_result%conf_high = min(delta + width, 1.0_dp)
      end if
   end function prop_test_real_vector

   pure function prop_test_int_vector(x, n, p, correct) result(test_result)
      !! Performs a test of equal or specified proportions for integer counts.
      integer, intent(in) :: x(:) !! Success counts for each sample.
      integer, intent(in) :: n(:) !! Trial counts conforming with `x`.
      real(dp), intent(in), optional :: p(:) !! Null probabilities conforming with `x`.
      logical, intent(in), optional :: correct !! Apply continuity correction for at most two samples.
      type(prop_test_result_t) :: test_result

      test_result = prop_test_real_vector(real(x, dp), real(n, dp), p, correct)
   end function prop_test_int_vector

   pure function cor_test_real(x, y, method) result(test_result)
      !! Performs a two-sided Pearson or asymptotic Spearman correlation test.
      real(dp), intent(in) :: x(:) !! First sample vector.
      real(dp), intent(in) :: y(:) !! Second sample vector conforming with `x`.
      character(len=*), intent(in), optional :: method !! `pearson` or `spearman`; defaults to Pearson.
      type(cor_test_result_t) :: test_result
      real(dp), allocatable :: ranks_x(:), ranks_y(:)
      real(dp) :: correlation, fisher_scale, fisher_z, t_statistic
      integer :: observation_count
      logical :: spearman

      observation_count = size(x)
      if (size(y) /= observation_count .or. observation_count < 3) return
      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) return
      spearman = .false.
      if (present(method)) then
         if (len(method) > 0) spearman = method(1:1) == "s" .or. method(1:1) == "S"
      end if
      if (spearman) then
         call r_average_ranks(x, ranks_x)
         call r_average_ranks(y, ranks_y)
         correlation = r_correlation(ranks_x, ranks_y)
         test_result%statistic = (real(observation_count, dp)**3 - real(observation_count, dp)) * &
                                 (1.0_dp - correlation) / 6.0_dp
         test_result%method = 2
      else
         correlation = r_correlation(x, y)
         test_result%method = 1
      end if
      test_result%estimate = correlation
      test_result%parameter = observation_count - 2
      if (ieee_is_nan(correlation)) return
      if (abs(correlation) >= 1.0_dp) return
      t_statistic = correlation * sqrt(real(test_result%parameter, dp) / &
                                       (1.0_dp - correlation**2))
      if (.not. spearman) test_result%statistic = t_statistic
      test_result%p_value = 2.0_dp * r_pt(abs(t_statistic), real(test_result%parameter, dp), &
                                         lower_tail=.false.)
      if (.not. spearman .and. observation_count > 3) then
         fisher_z = atanh(correlation)
         fisher_scale = r_qnorm(0.975_dp) / sqrt(real(observation_count - 3, dp))
         test_result%conf_low = tanh(fisher_z - fisher_scale)
         test_result%conf_high = tanh(fisher_z + fisher_scale)
      end if
   end function cor_test_real

   pure function cor_test_integer(x, y, method) result(test_result)
      !! Performs a Pearson or asymptotic Spearman test for integer samples.
      integer, intent(in) :: x(:) !! First sample vector.
      integer, intent(in) :: y(:) !! Second sample vector conforming with `x`.
      character(len=*), intent(in), optional :: method !! `pearson` or `spearman`; defaults to Pearson.
      type(cor_test_result_t) :: test_result

      test_result = cor_test_real(real(x, dp), real(y, dp), method)
   end function cor_test_integer

   pure function t_test_one(x, mu) result(test_result)
      !! Performs a two-sided one-sample Student t test and returns a 95 percent confidence interval.
      real(dp), intent(in) :: x(:) !! Sample observations.
      real(dp), intent(in), optional :: mu !! Mean under the null hypothesis; defaults to zero.
      type(t_test_result_t) :: test_result
      real(dp) :: critical_value, sample_mean, sample_variance
      integer :: sample_size

      sample_size = size(x)
      test_result%null_value = optval(mu, 0.0_dp)
      if (sample_size <= 1) return

      sample_mean = sum(x) / real(sample_size, dp)
      sample_variance = sum((x - sample_mean)**2) / real(sample_size - 1, dp)
      test_result%stderr = sqrt(sample_variance / real(sample_size, dp))
      test_result%estimate = sample_mean
      test_result%parameter = real(sample_size - 1, dp)
      test_result%method = 1
      if (test_result%stderr <= 0.0_dp) return

      test_result%statistic = (sample_mean - test_result%null_value) / test_result%stderr
      test_result%p_value = 2.0_dp * r_pt(abs(test_result%statistic), &
                                         test_result%parameter, lower_tail=.false.)
      critical_value = r_qt(0.975_dp, test_result%parameter)
      test_result%conf_low = sample_mean - critical_value * test_result%stderr
      test_result%conf_high = sample_mean + critical_value * test_result%stderr
   end function t_test_one

   pure function t_test_two(x, y, paired, var_equal) result(test_result)
      !! Performs a two-sided independent or paired two-sample t test with a 95 percent interval.
      real(dp), intent(in) :: x(:) !! First sample observations.
      real(dp), intent(in) :: y(:) !! Second sample observations.
      logical, intent(in), optional :: paired !! Use paired differences when true; defaults to false.
      logical, intent(in), optional :: var_equal !! Pool independent-sample variances when true.
      type(t_test_result_t) :: test_result
      real(dp), allocatable :: differences(:)
      real(dp) :: critical_value, denominator, mean_x, mean_y, n_x, n_y
      real(dp) :: pooled_variance, variance_x, variance_y

      if (optval(paired, .false.)) then
         if (size(x) /= size(y)) return
         differences = x - y
         test_result = t_test_one(differences)
         test_result%method = 4
         return
      end if

      if (size(x) <= 1 .or. size(y) <= 1) return
      n_x = real(size(x), dp)
      n_y = real(size(y), dp)
      mean_x = sum(x) / n_x
      mean_y = sum(y) / n_y
      variance_x = sum((x - mean_x)**2) / (n_x - 1.0_dp)
      variance_y = sum((y - mean_y)**2) / (n_y - 1.0_dp)
      test_result%estimate = mean_x
      test_result%estimate2 = mean_y

      if (optval(var_equal, .false.)) then
         pooled_variance = ((n_x - 1.0_dp) * variance_x + (n_y - 1.0_dp) * variance_y) / &
                           (n_x + n_y - 2.0_dp)
         test_result%stderr = sqrt(pooled_variance * (1.0_dp / n_x + 1.0_dp / n_y))
         test_result%parameter = n_x + n_y - 2.0_dp
         test_result%method = 3
      else
         test_result%stderr = sqrt(variance_x / n_x + variance_y / n_y)
         denominator = (variance_x / n_x)**2 / (n_x - 1.0_dp) + &
                       (variance_y / n_y)**2 / (n_y - 1.0_dp)
         if (denominator > 0.0_dp) then
            test_result%parameter = (variance_x / n_x + variance_y / n_y)**2 / denominator
         end if
         test_result%method = 2
      end if
      if (test_result%stderr <= 0.0_dp) return

      test_result%statistic = (mean_x - mean_y) / test_result%stderr
      test_result%p_value = 2.0_dp * r_pt(abs(test_result%statistic), &
                                         test_result%parameter, lower_tail=.false.)
      critical_value = r_qt(0.975_dp, test_result%parameter)
      test_result%conf_low = mean_x - mean_y - critical_value * test_result%stderr
      test_result%conf_high = mean_x - mean_y + critical_value * test_result%stderr
   end function t_test_two

   pure function t_test_p_value_one(x, mu) result(p_value)
      !! Returns only the two-sided p-value from a one-sample t test.
      real(dp), intent(in) :: x(:) !! Sample observations.
      real(dp), intent(in), optional :: mu !! Mean under the null hypothesis; defaults to zero.
      real(dp) :: p_value
      type(t_test_result_t) :: test_result

      test_result = t_test_one(x, mu)
      p_value = test_result%p_value
   end function t_test_p_value_one

   pure function t_test_p_value_two(x, y, paired, var_equal) result(p_value)
      !! Returns only the two-sided p-value from an independent or paired two-sample t test.
      real(dp), intent(in) :: x(:) !! First sample observations.
      real(dp), intent(in) :: y(:) !! Second sample observations.
      logical, intent(in), optional :: paired !! Use paired differences when true; defaults to false.
      logical, intent(in), optional :: var_equal !! Pool independent-sample variances when true.
      real(dp) :: p_value
      type(t_test_result_t) :: test_result

      test_result = t_test_two(x, y, paired, var_equal)
      p_value = test_result%p_value
   end function t_test_p_value_two

end module r_stats_tests
