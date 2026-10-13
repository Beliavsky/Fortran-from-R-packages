! SPDX-License-Identifier: GPL-2.0-or-later
program test_nortest
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan, ieee_is_finite
   use nortest
   use nortest_reference
   implicit none
   real(dp), allocatable :: x(:), changed(:)
   real(dp) :: nan, inf
   type(normality_result) :: results(8), baseline(8), altered(8), item
   integer :: sample, i, assertions

   assertions = 0
   do sample = 1, reference_count
      call make_sample(sample, x)
      results = evaluate(x)
      do i = 1, size(results)
         call check(results(i)%n == reference_n(sample), 'reference sample size')
         call check(results(i)%status == reference_status(i, sample), 'reference status')
         if (results(i)%status /= normality_success) then
            call check(ieee_is_nan(results(i)%statistic), 'error statistic is NaN')
            call check(ieee_is_nan(results(i)%p_value), 'error probability is NaN')
            cycle
         end if
         call close(results(i)%statistic, reference_statistic(i, sample), 5.0e-11_dp, 'R statistic', sample, i)
         call probability_close(results(i)%p_value, reference_p_value(i, sample), sample, i)
         call check(results(i)%n_classes == reference_classes(i, sample), 'R class count')
         call check(results(i)%df == reference_df(i, sample), 'R degrees of freedom')
      end do
      if (sample == 11) then
         call check(results(1)%p_value_at_floor, 'AD lower approximation limit')
         call check(results(2)%p_value_at_floor, 'CvM lower approximation limit')
         call close(results(1)%p_value, 3.7e-24_dp, 0.0_dp, 'AD floor', sample, 1)
         call close(results(2)%p_value, 7.37e-10_dp, 0.0_dp, 'CvM floor', sample, 2)
      end if
   end do

   call make_sample(1, x)
   baseline = evaluate(x)
   changed = x
   altered = evaluate(changed)
   call check(all(changed == x), 'input data unchanged')
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   changed = [nan, x(20:1:-1), nan]
   altered = evaluate(changed)
   call same_results(baseline, altered, 'NaN omission and reversed order')
   changed = 13.0_dp - 2.5_dp * x
   altered = evaluate(changed)
   call same_results(baseline, altered, 'negative affine transformation')
   changed = x * 1.0e250_dp
   altered = evaluate(changed)
   call same_results(baseline, altered, 'large finite scaling')
   changed = x * 1.0e-250_dp
   altered = evaluate(changed)
   call same_results(baseline, altered, 'small finite scaling')

   altered = evaluate([real(dp) ::])
   call check(all(altered%status == normality_invalid_size), 'empty sample rejected')
   altered = evaluate([nan, nan])
   call check(all(altered%status == normality_invalid_size), 'all missing sample rejected')
   altered = evaluate([1.0_dp])
   call check(all(altered%status == normality_invalid_size), 'singleton sample rejected')
   altered = evaluate(spread(3.0_dp, 1, 20))
   call check(all(altered%status == normality_zero_variance), 'constant sample rejected')
   altered = evaluate(spread(0.1_dp, 1, 20))
   call check(all(altered%status == normality_zero_variance), 'noninteger constant sample rejected')
   altered = evaluate(spread(0.0_dp, 1, 20))
   call check(all(altered%status == normality_zero_variance), 'zero sample rejected')
   altered = evaluate([x, inf])
   call check(all(altered%status == normality_nonfinite), 'infinity rejected')
   altered = evaluate([x, -inf])
   call check(all(altered%status == normality_nonfinite), 'negative infinity rejected')
   altered = evaluate(spread(huge(1.0_dp), 1, 20))
   call check(all(altered%status == normality_zero_variance), 'huge constant sample rejected')
   item = sf_test(spread(1.0_dp, 1, 5001))
   call check(item%status == normality_invalid_size, 'SF upper size bound')
   item = ad_test(x(:7))
   call check(item%status == normality_invalid_size, 'AD minimum size')
   item = cvm_test(x(:7))
   call check(item%status == normality_invalid_size, 'CvM minimum size')
   item = lillie_test(x(:4))
   call check(item%status == normality_invalid_size, 'Lilliefors minimum size')
   item = sf_test(x(:4))
   call check(item%status == normality_invalid_size, 'SF minimum size')
   item = pearson_test(x, n_classes = 0)
   call check(item%status == normality_invalid_classes, 'zero classes rejected')
   item = pearson_test(x, n_classes = -1)
   call check(item%status == normality_invalid_classes, 'negative classes rejected')
   item = pearson_test(x, n_classes = 2)
   call check(item%status == normality_invalid_classes, 'negative degrees of freedom rejected')
   item = pearson_test([-1.0_dp, 1.0_dp], n_classes = 1, adjust = .false.)
   call check(item%status == normality_success .and. item%df == 0, 'zero degrees of freedom accepted')
   call close(item%statistic, 0.0_dp, 0.0_dp, 'zero df statistic', 0, 0)
   call close(item%p_value, 1.0_dp, 0.0_dp, 'R zero df zero statistic p-value', 0, 0)
   item = pearson_test([-1.0_dp, 1.0_dp])
   call check(item%status == normality_success .and. item%df == 0, 'two observation default Pearson')
   call close(item%statistic, 1.0_dp, 1.0e-14_dp, 'two observation Pearson statistic', 0, 0)
   call close(item%p_value, 0.0_dp, 0.0_dp, 'R zero df positive statistic p-value', 0, 0)
   item = pearson_test(x, n_classes = huge(1))
   call check(item%status == normality_success, 'large class count without dense allocation')
   call check(ieee_is_finite(item%statistic) .and. ieee_is_finite(item%p_value), 'large class count finite result')
   print '(a, i0, a)', 'Passed ', assertions, ' nortest assertions (104 upstream R evaluations).'

contains

   pure function evaluate(data) result(tests)
      real(dp), intent(in) :: data(:) !! Sample passed to all five normality tests and three extra Pearson variants.
      type(normality_result) :: tests(8)

      tests = [ad_test(data), cvm_test(data), lillie_test(data), pearson_test(data), &
               pearson_test(data, adjust = .false.), pearson_test(data, n_classes = 7), &
               pearson_test(data, n_classes = 7, adjust = .false.), sf_test(data)]
   end function evaluate

   pure subroutine make_sample(index, data)
      integer, intent(in) :: index !! One-based fixture number matching test/reference/generate.R.
      real(dp), allocatable, intent(out) :: data(:) !! Deterministic sample corresponding to the selected R fixture.
      real(dp) :: nr, probability, angle
      integer :: j

      allocate(data(reference_n(index)))
      select case (index)
      case (1)
         data = [-1.6_dp, -1.1_dp, -0.9_dp, -0.8_dp, -0.6_dp, -0.5_dp, -0.4_dp, -0.3_dp, -0.2_dp, -0.1_dp, &
                 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.7_dp, 0.9_dp, 1.2_dp, 1.8_dp]
      case (2)
         data = [(sin(real(j, dp) * 1.7_dp) + cos(real(j, dp) * 0.31_dp), j = 1, size(data))]
      case (3)
         data = normal_quantiles
      case (4)
         data = normal_quantiles + 0.15_dp * normal_quantiles ** 2
      case (5)
         data = normal_quantiles + 0.25_dp * normal_quantiles ** 2
      case (6)
         data = normal_quantiles + 0.30_dp * normal_quantiles ** 2
      case (7, 8)
         data = [(real(j, dp), j = 1, size(data))]
      case (9, 11)
         data = 0.0_dp
         data(size(data)) = 1.0_dp
      case (10)
         data(:100) = -1.0_dp
         data(101:) = 1.0_dp
      case (12)
         nr = acos(-1.0_dp)
         do j = 1, size(data)
            probability = (real(mod(j * 37, 997), dp) + 0.5_dp) / 997.0_dp
            angle = 2.0_dp * nr * (real(mod(j * 67, 991), dp) + 0.5_dp) / 991.0_dp
            data(j) = sqrt(-2.0_dp * log(probability)) * cos(angle)
         end do
      case (13)
         data = [-2.0_dp, -0.7_dp, 0.1_dp, 0.8_dp, 1.7_dp]
      end select
   end subroutine make_sample

   subroutine same_results(expected, actual, label)
      type(normality_result), intent(in) :: expected(:) !! Baseline results before reordering, adding NaNs, or scaling.
      type(normality_result), intent(in) :: actual(:) !! Results of the transformed but statistically equivalent sample.
      character(len = *), intent(in) :: label !! Name of the invariance property being checked.
      integer :: j

      do j = 1, size(expected)
         call check(actual(j)%status == expected(j)%status, label)
         call check(actual(j)%n == expected(j)%n, label)
         call close(actual(j)%statistic, expected(j)%statistic, 5.0e-11_dp, label, 0, j)
         call close(actual(j)%p_value, expected(j)%p_value, 5.0e-10_dp, label, 0, j)
      end do
   end subroutine same_results

   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion that must hold for the test suite to pass.
      character(len = *), intent(in) :: label !! Diagnostic printed if the assertion fails.

      assertions = assertions + 1
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         error stop 1
      end if
   end subroutine check

   subroutine probability_close(actual, expected, sample_index, test_index)
      real(dp), intent(in) :: actual !! Computed p-value, including representable subnormal tail probabilities.
      real(dp), intent(in) :: expected !! Frozen upstream R p-value for the same sample and test variant.
      integer, intent(in) :: sample_index !! Fixture number used in the error diagnostic.
      integer, intent(in) :: test_index !! Test variant used in the error diagnostic.

      if (expected == 0.0_dp) then
         call close(actual, expected, 0.0_dp, 'R underflowed p-value', sample_index, test_index)
      else if (expected < 1.0e-6_dp) then
         call check(ieee_is_finite(actual) .and. actual > 0.0_dp, 'positive representable R tail retained')
         call close(log(actual / expected), 0.0_dp, 5.0e-8_dp, 'R log-ratio tail p-value', sample_index, test_index)
      else
         call close(actual, expected, 5.0e-10_dp, 'R p-value', sample_index, test_index)
      end if
   end subroutine probability_close

   subroutine close(actual, expected, tolerance, label, sample_index, test_index)
      real(dp), intent(in) :: actual !! Computed numerical result, which must be finite.
      real(dp), intent(in) :: expected !! Frozen R value or independently known identity.
      real(dp), intent(in) :: tolerance !! Absolute tolerance scaled by max(1, abs(expected)).
      character(len = *), intent(in) :: label !! Numerical quantity or property under test.
      integer, intent(in) :: sample_index !! Fixture number for diagnostics, or zero for an additional edge test.
      integer, intent(in) :: test_index !! Normality test variant for diagnostics, or zero for a standalone test.

      assertions = assertions + 1
      if (.not. ieee_is_finite(actual)) then
         print '(a)', 'Nonfinite result: ' // label
         error stop 1
      end if
      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         print '(a, 2i6)', 'FAIL: ' // label, sample_index, test_index
         print '(a, es25.17)', 'actual:   ', actual
         print '(a, es25.17)', 'expected: ', expected
         error stop 1
      end if
   end subroutine close

end program test_nortest
