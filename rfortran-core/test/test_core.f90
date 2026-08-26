! SPDX-License-Identifier: MIT
program test_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_mod, only : dp, r_autocorrelation, r_autocovariance, r_count_nonmissing, &
      r_cross_correlation, r_cross_covariance, r_difference, r_invalid_input, &
      r_average_ranks, r_dnorm, r_log_mean_exp, r_log_sum_exp, r_mean, r_ok, &
      r_pnorm, r_quantile_type7, r_sd, r_singular, r_variance, r_digamma, &
      r_log_beta, r_roll_mean_right, r_roll_mean_valid, r_trigamma
   implicit none

   real(dp) :: x(5), with_na(3), constant(4), matrix_data(8, 2), difference_input(4)
   real(dp), allocatable :: values(:)
   real(dp), allocatable :: matrix_values(:,:,:), differences(:), matrix_differences(:,:), ranks(:)
   integer :: i, status

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   call check_close(r_mean(x), 3.0_dp, 1.0e-14_dp, 'mean')
   call check_close(r_variance(x), 2.5_dp, 1.0e-14_dp, 'sample variance')
   call check_close(r_variance(x, ddof=0), 2.0_dp, 1.0e-14_dp, 'population variance')
   call check_close(r_sd(x), sqrt(2.5_dp), 1.0e-14_dp, 'sample standard deviation')

   with_na = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp]
   call check(ieee_is_nan(r_mean(with_na)), 'NaN propagation')
   call check_close(r_mean(with_na, na_rm=.true.), 2.0_dp, 1.0e-14_dp, 'NaN removal')
   call check(r_count_nonmissing(with_na) == 2, 'nonmissing count')
   call check(ieee_is_nan(r_variance([1.0_dp])), 'insufficient variance data')

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

   constant = 1.0_dp
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
   call r_roll_mean_valid(x, 3, values, status)
   call check(status == r_ok .and. all(abs(values - [2.0_dp, 3.0_dp, 4.0_dp]) < 1.0e-14_dp), &
      'valid rolling mean')
   call r_roll_mean_right(x, 3, values, status)
   call check(status == r_ok .and. all(ieee_is_nan(values(:2))) .and. &
      all(abs(values(3:) - [2.0_dp, 3.0_dp, 4.0_dp]) < 1.0e-14_dp), 'right-aligned rolling mean')
   call r_roll_mean_valid([1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp, 4.0_dp], 2, values, status)
   call check(ieee_is_nan(values(1)) .and. ieee_is_nan(values(2)) .and. &
      abs(values(3) - 3.5_dp) < 1.0e-14_dp, 'rolling mean recovers after NaN leaves window')

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
