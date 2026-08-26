! SPDX-License-Identifier: MIT
program test_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_mod, only : dp, r_autocorrelation, r_autocovariance, r_count_nonmissing, &
      r_cross_correlation, r_cross_covariance, r_difference, r_invalid_input, &
      r_average_ranks, r_correlation, r_covariance, r_dnorm, r_log_mean_exp, r_log_sum_exp, r_mean, r_ok, &
      r_pnorm, r_qnorm, r_quantile_type7, r_sd, r_singular, r_variance, r_digamma, &
      r_log_beta, r_log_choose, r_log_factorial, r_roll_correlation_right, r_roll_covariance_valid, &
      r_roll_max_valid, r_roll_mean_right, r_roll_mean_valid, r_roll_min_right, r_roll_min_valid, &
      r_roll_sd_right, r_roll_sum_valid, r_roll_variance_valid, r_trigamma
   use r_mod, only : r_weighted_correlation, r_weighted_covariance, r_weighted_mean
   use r_mod, only : r_qrule_hf7, r_qrule_school
   use r_mod, only : r_weighted_quantile_ecdf, r_weighted_quantile_frequency_type7
   use r_mod, only : r_weighted_quantile_isotone, r_weighted_quantile_linear_cdf
   use r_mod, only : r_weighted_quantile_survey, r_weighted_sd, r_weighted_variance
   use r_mod, only : r_mad, r_median, r_order
   use r_mod, only : r_df, r_dchisq, r_dt, r_pf, r_pchisq, r_pt
   use r_mod, only : r_qf, r_qchisq, r_qt
   use r_mod, only : r_regularized_beta, r_regularized_gamma_p, r_regularized_gamma_q
   use r_transforms, only : r_expm1, r_log1mexp, r_log1p, r_log1pexp, r_logistic, r_logit
   implicit none

   real(dp) :: x(5), with_na(3), constant(4), matrix_data(8, 2), difference_input(4)
   real(dp), allocatable :: values(:)
   real(dp), allocatable :: matrix_values(:,:,:), differences(:), matrix_differences(:,:), ranks(:)
   integer, allocatable :: order(:)
   integer :: i, status

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   call check_close(r_mean(x), 3.0_dp, 1.0e-14_dp, 'mean')
   call check_close(r_variance(x), 2.5_dp, 1.0e-14_dp, 'sample variance')
   call check_close(r_variance(x, ddof=0), 2.0_dp, 1.0e-14_dp, 'population variance')
   call check_close(r_sd(x), sqrt(2.5_dp), 1.0e-14_dp, 'sample standard deviation')

   with_na = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp]
   constant = 1.0_dp
   call check(ieee_is_nan(r_mean(with_na)), 'NaN propagation')
   call check_close(r_mean(with_na, na_rm=.true.), 2.0_dp, 1.0e-14_dp, 'NaN removal')
   call check(r_count_nonmissing(with_na) == 2, 'nonmissing count')
   call check(ieee_is_nan(r_variance([1.0_dp])), 'insufficient variance data')
   call check_close(r_covariance([1.0_dp, 2.0_dp, 3.0_dp], [2.0_dp, 4.0_dp, 6.0_dp]), &
      2.0_dp, 1.0e-14_dp, 'sample covariance')
   call check_close(r_correlation([1.0_dp, 2.0_dp, 3.0_dp], [2.0_dp, 4.0_dp, 6.0_dp]), &
      1.0_dp, 1.0e-14_dp, 'correlation')
   call check(ieee_is_nan(r_correlation(constant, constant)), 'constant correlation')
   call check_close(r_covariance(with_na, [2.0_dp, 9.0_dp, 6.0_dp], na_rm=.true.), &
      4.0_dp, 1.0e-14_dp, 'pairwise covariance NaN removal')
   call check_close(r_covariance([1.0_dp, ieee_value(0.0_dp, ieee_positive_inf), 3.0_dp], &
      [2.0_dp, 9.0_dp, 6.0_dp], finite_only=.true.), 4.0_dp, 1.0e-14_dp, 'pairwise finite filtering')
   call check_close(r_weighted_mean([1.0_dp, 2.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp]), &
      17.0_dp/6.0_dp, 1.0e-14_dp, 'weighted mean')
   call check_close(r_weighted_variance([1.0_dp, 2.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp]), &
      53.0_dp/36.0_dp, 1.0e-14_dp, 'weighted population variance')
   call check_close(r_weighted_variance([1.0_dp, 2.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp], unbiased=.true.), &
      53.0_dp/22.0_dp, 1.0e-14_dp, 'weighted unbiased variance')
   call check_close(r_weighted_sd([1.0_dp, 2.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 3.0_dp]), &
      sqrt(53.0_dp/36.0_dp), 1.0e-14_dp, 'weighted standard deviation')
   call check_close(r_weighted_covariance([1.0_dp, 2.0_dp, 4.0_dp], [2.0_dp, 1.0_dp, 8.0_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp], unbiased=.true.), 68.0_dp/11.0_dp, 1.0e-14_dp, 'weighted covariance')
   call check_close(r_weighted_correlation([1.0_dp, 2.0_dp, 4.0_dp], [2.0_dp, 1.0_dp, 8.0_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp]), (34.0_dp/9.0_dp)/sqrt((53.0_dp/36.0_dp)*(101.0_dp/9.0_dp)), &
      1.0e-14_dp, 'weighted correlation')
   call check_close(r_weighted_mean(with_na, [1.0_dp, 0.0_dp, 1.0_dp]), 2.0_dp, 1.0e-14_dp, &
      'zero weight excludes NaN')
   call check(ieee_is_nan(r_weighted_mean(x, [-1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])), &
      'negative weight rejection')
   call check_close(r_weighted_quantile_ecdf([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.5_dp), 2.0_dp, 1.0e-14_dp, 'weighted ECDF median')
   call check_close(r_weighted_quantile_ecdf([1.0_dp, 2.0_dp, 3.0_dp], &
      [0.0_dp, 1.0_dp, 3.0_dp], 0.5_dp), 3.0_dp, 1.0e-14_dp, 'weighted ECDF threshold')
   call check_close(r_weighted_quantile_ecdf(with_na, [1.0_dp, 0.0_dp, 1.0_dp], 0.5_dp), &
      1.0_dp, 1.0e-14_dp, 'weighted quantile excludes zero-weight NaN')
   call check(ieee_is_nan(r_weighted_quantile_ecdf(x, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      -0.1_dp)), 'weighted quantile probability validation')
   call check_close(r_weighted_quantile_linear_cdf([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.625_dp), 2.5_dp, 1.0e-14_dp, &
      'weighted cumulative-probability interpolation')
   call check_close(r_weighted_quantile_frequency_type7([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.5_dp), 3.0_dp, 1.0e-14_dp, &
      'frequency-weighted type-7 convention')
   call check_close(r_weighted_quantile_isotone([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.5_dp), 2.5_dp, 1.0e-14_dp, &
      'isotone exact-boundary interpolation')
   call check_close(r_weighted_quantile_survey([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.5_dp, r_qrule_school), &
      2.5_dp, 1.0e-14_dp, 'survey school rule')
   call check_close(r_weighted_quantile_survey([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.25_dp, r_qrule_hf7), &
      1.75_dp, 1.0e-14_dp, 'survey Hyndman-Fan type 7 rule')
   call check_close(r_median([5.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]), &
      3.0_dp, 1.0e-14_dp, 'median')
   call check_close(r_mad([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]), &
      1.4826_dp, 1.0e-14_dp, 'median absolute deviation')
   call r_order([2.0_dp, 1.0_dp, 2.0_dp], order)
   call check(all(order == [2, 1, 3]), 'stable ordering')

   call r_autocovariance(x, values, lag_max=2, status=status)
   call check(status == r_ok, 'autocovariance status')
   call check(lbound(values, 1) == 0 .and. ubound(values, 1) == 2, 'lag bounds')
   call check_close(values(0), 2.0_dp, 1.0e-14_dp, 'autocovariance lag zero')
   call check_close(values(1), 0.8_dp, 1.0e-14_dp, 'autocovariance lag one')
   call check_close(values(2), -0.2_dp, 1.0e-14_dp, 'autocovariance lag two')

   call r_autocorrelation(x, values, lag_max=2, status=status)
   call check(status == r_ok, 'autocorrelation status')
   call check_close(values(0), 1.0_dp, 1.0e-14_dp, 'autocorrelation lag zero')
   call check_close(values(1), 0.4_dp, 1.0e-14_dp, 'autocorrelation lag one')
   call check_close(values(2), -0.1_dp, 1.0e-14_dp, 'autocorrelation lag two')

   call r_autocorrelation(constant, values, lag_max=1, status=status)
   call check(status == r_singular, 'constant autocorrelation status')
   call r_autocovariance(with_na, values, lag_max=1, status=status)
   call check(status == r_invalid_input, 'NaN autocovariance status')

   do i = 1, 8
      matrix_data(i, 1) = real(i, dp)
      matrix_data(i, 2) = 2.0_dp*real(i, dp) + 1.0_dp
   end do
   call r_cross_covariance(matrix_data, matrix_values, lag_max=2, status=status)
   call check(status == r_ok, 'cross covariance status')
   call check_close(matrix_values(0, 1, 1), 5.25_dp, 1.0e-14_dp, 'cross covariance lag zero')
   call r_cross_correlation(matrix_data, matrix_values, lag_max=2, status=status)
   call check(status == r_ok, 'cross correlation status')
   call check_close(matrix_values(0, 1, 2), 1.0_dp, 1.0e-14_dp, 'cross correlation lag zero')

   difference_input = [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp]
   call r_difference(difference_input, differences, status=status)
   call check(status == r_ok .and. size(differences) == 3, 'vector difference status and shape')
   call check(all(abs(differences - [1.0_dp, 2.0_dp, 3.0_dp]) < 1.0e-14_dp), 'vector difference values')
   call r_difference(difference_input, differences, differences=2, status=status)
   call check(all(abs(differences - [1.0_dp, 1.0_dp]) < 1.0e-14_dp), 'second differences')
   call r_difference(matrix_data, matrix_differences, status=status)
   call check(status == r_ok .and. all(shape(matrix_differences) == [7, 2]), 'matrix difference shape')
   call check(all(abs(matrix_differences(:, 2) - 2.0_dp) < 1.0e-14_dp), 'matrix difference values')

   call check_close(r_dnorm(0.0_dp), 1.0_dp/sqrt(2.0_dp*acos(-1.0_dp)), 1.0e-14_dp, 'normal density')
   call check_close(r_pnorm(0.0_dp), 0.5_dp, 1.0e-14_dp, 'normal CDF')
   call check_close(r_pnorm(10.0_dp, lower_tail=.false., log_probability=.true.), &
      -53.23128515051247_dp, 2.0e-13_dp, 'normal log upper tail')
   call check_close(r_qnorm(0.975_dp), 1.959963984540054_dp, 2.0e-14_dp, 'normal quantile')
   call check_close(r_qnorm(log(0.025_dp), log_probability=.true.), &
      -1.959963984540054_dp, 2.0e-14_dp, 'normal log-probability quantile')
   call check_close(r_qnorm(log(0.025_dp), lower_tail=.false., log_probability=.true.), &
      1.959963984540054_dp, 2.0e-14_dp, 'normal log upper-tail quantile')
   call check_close(r_qnorm(0.975_dp, location=2.0_dp, scale=3.0_dp), &
      2.0_dp + 3.0_dp*1.959963984540054_dp, 2.0e-14_dp, 'scaled normal quantile')
   call check_close(r_pnorm(r_qnorm(-1000.0_dp, log_probability=.true.), &
      log_probability=.true.), -1000.0_dp, 2.0e-13_dp, 'extreme normal log-tail inversion')
   call check(ieee_is_nan(r_qnorm(-0.1_dp)), 'normal quantile invalid probability')
   call check_close(r_quantile_type7([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], 0.25_dp), &
      1.75_dp, 1.0e-14_dp, 'type-7 quantile')
   call r_average_ranks([10.0_dp, 20.0_dp, 20.0_dp, 30.0_dp], ranks, status=status)
   call check(status == r_ok, 'average rank status')
   call check(all(abs(ranks - [1.0_dp, 2.5_dp, 2.5_dp, 4.0_dp]) < 1.0e-14_dp), 'average ranks')
   call check_close(r_log_sum_exp([1000.0_dp, 1000.0_dp]), 1000.0_dp + log(2.0_dp), &
      1.0e-14_dp, 'log-sum-exp')
   call check_close(r_log_mean_exp([1000.0_dp, 1000.0_dp]), 1000.0_dp, 1.0e-14_dp, 'log-mean-exp')
   call check_close(r_digamma(1.0_dp), -0.5772156649015328606_dp, 2.0e-12_dp, 'digamma')
   call check_close(r_trigamma(1.0_dp), acos(-1.0_dp)**2/6.0_dp, 2.0e-12_dp, 'trigamma')
   call check_close(r_log_beta(2.0_dp, 3.0_dp), log(1.0_dp/12.0_dp), 1.0e-14_dp, 'log beta')
   call check_close(r_regularized_gamma_p(2.0_dp, 3.0_dp), &
      1.0_dp - 4.0_dp*exp(-3.0_dp), 2.0e-14_dp, 'regularized lower gamma')
   call check_close(r_regularized_gamma_q(2.0_dp, 3.0_dp), &
      4.0_dp*exp(-3.0_dp), 2.0e-14_dp, 'regularized upper gamma')
   call check_close(r_regularized_beta(0.5_dp, 2.0_dp, 3.0_dp), &
      0.6875_dp, 2.0e-14_dp, 'regularized beta')
   call check_close(r_dt(0.0_dp,1.0_dp),1.0_dp/acos(-1.0_dp),2.0e-14_dp,'Student-t density')
   call check_close(r_pt(1.0_dp,1.0_dp),0.75_dp,2.0e-14_dp,'Student-t CDF')
   call check_close(exp(r_pt(1.0_dp,1.0_dp,lower_tail=.false.,log_probability=.true.)), &
      0.25_dp,2.0e-14_dp,'Student-t log survival')
   call check_close(r_qt(0.75_dp,1.0_dp),1.0_dp,2.0e-13_dp,'Student-t quantile')
   call check_close(r_pt(r_qt(0.95_dp,5.0_dp),5.0_dp),0.95_dp,2.0e-13_dp,'Student-t round trip')
   call check_close(r_dchisq(2.0_dp,2.0_dp),0.5_dp*exp(-1.0_dp),2.0e-14_dp,'chi-square density')
   call check_close(r_pchisq(2.0_dp,2.0_dp),1.0_dp-exp(-1.0_dp),2.0e-14_dp,'chi-square CDF')
   call check_close(r_qchisq(1.0_dp-exp(-1.0_dp),2.0_dp),2.0_dp,2.0e-13_dp,'chi-square quantile')
   call check_close(r_df(1.0_dp,2.0_dp,2.0_dp),0.25_dp,2.0e-14_dp,'F density')
   call check_close(r_pf(1.0_dp,2.0_dp,2.0_dp),0.5_dp,2.0e-14_dp,'F CDF')
   call check_close(r_qf(0.75_dp,2.0_dp,2.0_dp),3.0_dp,2.0e-13_dp,'F quantile')
   call check_close(r_pf(r_qf(0.95_dp,5.0_dp,9.0_dp),5.0_dp,9.0_dp),0.95_dp,2.0e-13_dp,'F round trip')
   call check(ieee_is_nan(r_pt(0.0_dp,-1.0_dp)),'Student-t invalid degrees of freedom')
   call check_close(r_log_factorial(5), log(120.0_dp), 1.0e-14_dp, 'log factorial')
   call check_close(r_log_choose(1000000, 2), log(499999500000.0_dp), 1.0e-14_dp, 'stable log choose')
   call check_close(r_log_choose(20, 3), r_log_choose(20, 17), 1.0e-14_dp, 'log choose symmetry')
   call check(ieee_is_nan(r_log_factorial(-1)), 'log factorial invalid domain')
   call r_roll_mean_valid(x, 3, values, status)
   call check(status == r_ok .and. all(abs(values - [2.0_dp, 3.0_dp, 4.0_dp]) < 1.0e-14_dp), &
      'valid rolling mean')
   call r_roll_mean_right(x, 3, values, status)
   call check(status == r_ok .and. all(ieee_is_nan(values(:2))) .and. &
      all(abs(values(3:) - [2.0_dp, 3.0_dp, 4.0_dp]) < 1.0e-14_dp), 'right-aligned rolling mean')
   call r_roll_mean_valid([1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp, 4.0_dp], 2, values, status)
   call check(ieee_is_nan(values(1)) .and. ieee_is_nan(values(2)) .and. &
      abs(values(3) - 3.5_dp) < 1.0e-14_dp, 'rolling mean recovers after NaN leaves window')
   call r_roll_sum_valid(x, 3, values, status=status)
   call check(status == r_ok .and. all(abs(values - [6.0_dp, 9.0_dp, 12.0_dp]) < 1.0e-14_dp), 'rolling sum')
   call r_roll_min_right(x, 3, values, status=status)
   call check(status == r_ok .and. all(ieee_is_nan(values(:2))) .and. &
      all(abs(values(3:) - [1.0_dp, 2.0_dp, 3.0_dp]) < 1.0e-14_dp), 'right-aligned rolling minimum')
   call r_roll_max_valid(x, 3, values, status=status)
   call check(status == r_ok .and. all(abs(values - [3.0_dp, 4.0_dp, 5.0_dp]) < 1.0e-14_dp), 'rolling maximum')
   call r_roll_sum_valid(with_na, 2, values, na_rm=.true., status=status)
   call check(all(abs(values - [1.0_dp, 3.0_dp]) < 1.0e-14_dp), 'rolling sum NaN removal')
   call r_roll_min_valid(with_na(2:2), 1, values, na_rm=.true., status=status)
   call check(.not. ieee_is_finite(values(1)) .and. values(1) > 0.0_dp, 'empty rolling minimum identity')
   call r_roll_max_valid(with_na(2:2), 1, values, na_rm=.true., status=status)
   call check(.not. ieee_is_finite(values(1)) .and. values(1) < 0.0_dp, 'empty rolling maximum identity')
   call r_roll_sum_valid([1.0_dp, ieee_value(0.0_dp, ieee_positive_inf), 3.0_dp], 3, &
      values, finite_only=.true., status=status)
   call check_close(values(1), 4.0_dp, 1.0e-14_dp, 'rolling sum finite filtering')
   call check_close(r_log1p(1.0e-12_dp), log(1.0_dp + 1.0e-12_dp), 1.0e-12_dp, 'log1p')
   call check_close(r_expm1(1.0e-12_dp), 1.0e-12_dp, 1.0e-14_dp, 'expm1')
   call check_close(r_log1mexp(-1.0_dp), log(1.0_dp - exp(-1.0_dp)), 1.0e-14_dp, 'log1mexp')
   call check_close(r_log1pexp(1000.0_dp), 1000.0_dp, 1.0e-14_dp, 'log1pexp upper tail')
   call check(ieee_is_nan(r_log1p(-2.0_dp)), 'log1p invalid domain')
   call check(ieee_is_nan(r_log1mexp(1.0_dp)), 'log1mexp invalid domain')
   call check_close(r_logistic(r_logit(0.25_dp)), 0.25_dp, 1.0e-14_dp, 'logistic/logit inverse')
   call check(ieee_is_nan(r_logit(-0.1_dp)), 'logit invalid domain')
   call r_roll_variance_valid(x, 3, values, status=status)
   call check(status == r_ok .and. all(abs(values - 1.0_dp) < 1.0e-14_dp), 'rolling sample variance')
   call r_roll_sd_right(x, 3, values, ddof=0, status=status)
   call check(status == r_ok .and. all(ieee_is_nan(values(:2))) .and. &
      all(abs(values(3:) - sqrt(2.0_dp/3.0_dp)) < 1.0e-14_dp), 'rolling population sd')
   call r_roll_variance_valid([1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp], 3, &
      values, na_rm=.true., status=status)
   call check_close(values(1), 2.0_dp, 1.0e-14_dp, 'rolling variance na_rm')
   call r_roll_covariance_valid(x, 2.0_dp*x + 1.0_dp, 3, values, status=status)
   call check(status == r_ok .and. all(abs(values - 2.0_dp) < 1.0e-14_dp), 'rolling sample covariance')
   call r_roll_correlation_right(x, 2.0_dp*x + 1.0_dp, 3, values, status=status)
   call check(status == r_ok .and. all(ieee_is_nan(values(:2))) .and. &
      all(abs(values(3:) - 1.0_dp) < 1.0e-14_dp), 'right-aligned rolling correlation')
   call r_roll_covariance_valid(x, x(:4), 3, values, status=status)
   call check(status == r_invalid_input .and. size(values) == 0, 'paired rolling shape validation')

   print '(a)', 'rfortran-core tests: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label

      call check(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected)), label)
   end subroutine check_close

end program test_core
