program test_quantile_types
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf
   use r_stats, only: dp, quantile, iqr
   use stats_test_assertions, only: assert_close, assert_vector_close, assert_true
   implicit none
   real(dp), parameter :: x(8) = [8.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 16.0_dp, 4.0_dp, 32.0_dp, 7.0_dp]
   real(dp), parameter :: probabilities(9) = [0.0_dp, 0.0625_dp, 0.125_dp, 0.25_dp, 0.3_dp, &
      0.5_dp, 0.75_dp, 0.9375_dp, 1.0_dp]
   real(dp) :: expected(9, 9), actual(9), nan, infinity
   integer :: kind

   expected(:, 1) = [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 32.0_dp, 32.0_dp]
   expected(:, 2) = [1.0_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.0_dp, 5.5_dp, 12.0_dp, 32.0_dp, 32.0_dp]
   expected(:, 3) = expected(:, 1)
   expected(:, 4) = [1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 24.0_dp, 32.0_dp]
   expected(:, 5) = expected(:, 2)
   expected(:, 6) = [1.0_dp, 1.0_dp, 1.125_dp, 2.0_dp, 2.0_dp, 5.5_dp, 14.0_dp, 32.0_dp, 32.0_dp]
   expected(:, 7) = [1.0_dp, 1.4375_dp, 1.875_dp, 2.0_dp, 2.2_dp, 5.5_dp, 10.0_dp, 25.0_dp, 32.0_dp]
   expected(:, 8) = [1.0_dp, 1.0_dp, 1.375_dp, 2.0_dp, 2.0_dp, 5.5_dp, &
      12.666666666666657_dp, 32.0_dp, 32.0_dp]
   expected(:, 9) = [1.0_dp, 1.0_dp, 1.40625_dp, 2.0_dp, 2.0_dp, 5.5_dp, 12.5_dp, 32.0_dp, 32.0_dp]
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   infinity = ieee_value(0.0_dp, ieee_positive_inf)
   do kind = 1, 9
      actual = quantile(x, probabilities, type=kind)
      call assert_true(all(ieee_is_finite(actual)), 'finite quantiles')
      call assert_vector_close(actual, expected(:, kind), 'R quantile types', 2.0e-12_dp)
      call assert_close(iqr(x, type=kind), expected(7, kind) - expected(4, kind), 'R IQR types', 2.0e-12_dp)
      call assert_close(quantile([x, nan], 0.75_dp, na_rm=.true., type=kind), &
         expected(7, kind), 'NaN removal', 2.0e-12_dp)
      call assert_true(ieee_is_nan(quantile([x, nan], 0.5_dp, type=kind)), 'NaN propagation')
      call assert_true(ieee_is_nan(quantile(x(:0), 0.5_dp, type=kind)), 'empty sample')
      call assert_true(ieee_is_nan(quantile(x, nan, type=kind)), 'NaN probability')
      call assert_close(quantile([3.0_dp], 0.4_dp, type=kind), 3.0_dp, 'singleton')
      call assert_true(quantile([infinity, infinity], 0.5_dp, type=kind) == infinity, 'equal infinities')
      call assert_true(quantile([1.0_dp, infinity], 0.0_dp, type=kind) == 1.0_dp, 'finite endpoint')
      call assert_true(quantile([1.0_dp, infinity], 1.0_dp, type=kind) == infinity, 'infinite endpoint')
   end do
   call assert_vector_close(quantile([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [0.375_dp, 0.625_dp, 0.875_dp], type=3), [2.0_dp, 2.0_dp, 4.0_dp], 'nearest even order statistic')
   print *, 'test_quantile_types: PASS'
end program test_quantile_types
