! SPDX-License-Identifier: MIT
program test_market
  use rtl, only: dp, beta_result, trade_stats_result, strategy_result
  use rtl, only: compute_returns, roll_adjust_mask, prompt_beta, trade_stats, moving_average_strategy
  implicit none

  real(dp) :: prices(5, 2), returns_matrix(6, 3), open_price(8), close_price(8)
  real(dp), allocatable :: ret(:, :)
  logical :: keep(5)
  type(beta_result) :: betas
  type(trade_stats_result) :: stats
  type(strategy_result) :: strategy

  prices(:, 1) = [10.0_dp, 11.0_dp, 13.0_dp, 12.0_dp, 15.0_dp]
  prices(:, 2) = [20.0_dp, 22.0_dp, 22.0_dp, 24.0_dp, 30.0_dp]
  ret = compute_returns(prices, 1, "abs")
  call assert_close(ret(1, 1), 1.0_dp, 0.0_dp)
  call assert_close(ret(4, 2), 6.0_dp, 0.0_dp)
  ret = compute_returns(prices, 1, "rel")
  call assert_close(ret(1, 1), 0.1_dp, 1.0e-14_dp)

  keep = roll_adjust_mask([1, 2, 3, 4, 5], [2, 4])
  call assert_true(all(keep .eqv. [.true., .true., .false., .true., .false.]))

  returns_matrix(:, 1) = [-0.03_dp, -0.02_dp, -0.01_dp, 0.01_dp, 0.02_dp, 0.03_dp]
  returns_matrix(:, 2) = 2.0_dp * returns_matrix(:, 1)
  returns_matrix(:, 3) = -returns_matrix(:, 1)
  betas = prompt_beta(returns_matrix)
  call assert_close(betas%all(1), 1.0_dp, 1.0e-14_dp)
  call assert_close(betas%all(2), 2.0_dp, 1.0e-14_dp)
  call assert_close(betas%bull(2), 2.0_dp, 1.0e-14_dp)
  call assert_close(betas%bear(3), -1.0_dp, 1.0e-14_dp)

  stats = trade_stats([0.01_dp, -0.02_dp, 0.03_dp, 0.0_dp], 0.0_dp, 4)
  call assert_close(stats%cumulative_return, product([1.01_dp, 0.98_dp, 1.03_dp, 1.0_dp]) - 1.0_dp, 1.0e-14_dp)
  call assert_close(stats%fraction_winning, 2.0_dp / 3.0_dp, 1.0e-14_dp)
  call assert_close(stats%fraction_in_market, 0.75_dp, 1.0e-14_dp)
  call assert_true(stats%maximum_drawdown <= 0.0_dp)

  open_price = [10.0_dp, 10.2_dp, 10.4_dp, 10.6_dp, 10.8_dp, 11.0_dp, 10.8_dp, 10.6_dp]
  close_price = [10.1_dp, 10.3_dp, 10.5_dp, 10.7_dp, 10.9_dp, 10.7_dp, 10.5_dp, 10.4_dp]
  strategy = moving_average_strategy(open_price, close_price, 2, 3)
  call assert_true(strategy%status%ok)
  call assert_true(size(strategy%signal) == 8)
  call assert_close(strategy%cumulative_equity(8), product(1.0_dp + strategy%strategy_return), 1.0e-14_dp)

  print '(a)', 'test_market: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print '(a,3es24.15)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_market
