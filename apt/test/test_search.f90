! SPDX-License-Identifier: GPL-2.0-or-later
program test_search
  use apt, only : dp, apt_tar, apt_mtar, ci_tar_lag_result, &
    ci_tar_threshold_result, ci_tar_lag, ci_tar_threshold
  use test_support, only : assert_true, assert_close, generate_prices
  implicit none
  real(dp) :: x(80), y(80)
  type(ci_tar_lag_result) :: lag_result
  type(ci_tar_threshold_result) :: threshold_result
  integer :: i
  call generate_prices(x, y)

  call ci_tar_lag(y, x, lag_result, apt_tar, 5, 0.0_dp, .true.)
  call assert_true(lag_result%status == 0, 'adjusted lag search status')
  do i = 2, size(lag_result%cointegration_observations)
    call assert_true(lag_result%cointegration_observations(i) == &
      lag_result%cointegration_observations(1), 'adjusted windows are equal')
  end do
  call assert_true(lag_result%best_lag_aic >= 0 .and. &
    lag_result%best_lag_aic <= 5, 'AIC lag range')
  call assert_close(lag_result%best_aic, minval(lag_result%aic), 1e-12_dp, 'minimum AIC')

  call ci_tar_threshold(y, x, threshold_result, apt_mtar, 1, 0.15_dp)
  call assert_true(threshold_result%status == 0, 'threshold search status')
  call assert_close(threshold_result%minimum_sse, minval(threshold_result%path_sse), &
    1e-12_dp, 'threshold minimum SSE')
  call assert_true(threshold_result%lower_index < threshold_result%upper_index, &
    'trimmed threshold range')
  call assert_true(any(abs(threshold_result%path_threshold - &
    threshold_result%threshold) < 1e-14_dp), 'selected threshold lies on path')
  print '(a)', 'test_search: PASS'
end program test_search
