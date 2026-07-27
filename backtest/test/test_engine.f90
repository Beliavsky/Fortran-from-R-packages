! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

program test_engine
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use backtest_kinds, only : dp
  use backtest_math, only : worst_drawdown
  use backtest_engine, only : backtest_config, backtest_result, &
                              run_numeric_backtest, run_grouped_backtest
  use backtest_summary, only : spread_statistics, sharpe_ratios, total_counts, &
                               marginal_counts, cumulative_bucket_returns
  implicit none

  real(dp) :: signals(12, 1), returns(12, 1), returns2(12, 2)
  integer :: period(12), ids(12), status
  integer :: grouped(12, 2)
  type(backtest_config) :: config
  type(backtest_result) :: result, grouped_result, by_result, universe_result
  real(dp), allocatable :: spread(:), ci_low(:), ci_high(:), sharpe(:), cumulative(:, :)
  real(dp) :: drawdown, return_loss
  integer :: draw_start, draw_end
  integer, allocatable :: totals(:, :), margins(:, :)
  logical :: universe(12)

  signals(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                   4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, &
                   1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]
  returns(:, 1) = [0.01_dp, 0.02_dp, 0.04_dp, 0.05_dp, &
                   0.06_dp, 0.04_dp, 0.01_dp, 0.00_dp, &
                   0.00_dp, 0.05_dp, 0.01_dp, 0.04_dp]
  period = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
  ids = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]

  config%n_buckets = 2
  config%by_period = .true.
  config%natural = .true.
  config%overlaps = 1
  call run_numeric_backtest(signals, returns, config, result, status, &
                            period=period, ids=ids)
  call assert_true(status == 0, "numeric backtest status")
  call assert_true(result%n_by == 3 .and. result%n_buckets == 2, "result dimensions")
  call assert_close(result%means(1, 1, 1, 1), 0.015_dp, 1.0e-12_dp, "period 1 low")
  call assert_close(result%means(1, 1, 1, 2), 0.045_dp, 1.0e-12_dp, "period 1 high")
  call assert_close(result%means(1, 1, 2, 1), 0.005_dp, 1.0e-12_dp, "period 2 low")
  call assert_close(result%means(1, 1, 2, 2), 0.05_dp, 1.0e-12_dp, "period 2 high")
  call assert_close(result%means(1, 1, 3, 1), 0.005_dp, 1.0e-12_dp, "period 3 low")
  call assert_close(result%means(1, 1, 3, 2), 0.045_dp, 1.0e-12_dp, "period 3 high")
  call assert_close(result%turnover(2, 1), 1.0_dp, 1.0e-12_dp, "complete portfolio replacement")
  call assert_close(result%turnover(3, 1), 0.5_dp, 1.0e-12_dp, "partial portfolio replacement")

  call spread_statistics(result, 1, 1, spread, ci_low, ci_high)
  call assert_array_close(spread, [0.03_dp, 0.045_dp, 0.04_dp], 1.0e-12_dp, &
                          "spread series")
  call assert_true(all(ieee_is_finite(ci_low)) .and. all(ieee_is_finite(ci_high)), &
                   "finite confidence bands")

  call sharpe_ratios(result, 1, sharpe)
  call assert_close(sharpe(1), mean3(spread) / sd3(spread), 1.0e-12_dp, "raw Sharpe ratio")

  call total_counts(result, 1, .true., totals)
  call assert_true(all(totals(:, 1) == 4), "low-high total counts")
  call marginal_counts(result, 1, 1, margins)
  call assert_true(margins(4, 3) == 12, "marginal grand total")

  call cumulative_bucket_returns(result, 1, 1, cumulative)
  call assert_close(cumulative(3, 1), 0.025_dp, 1.0e-12_dp, "cumulative low return")
  call assert_close(cumulative(3, 2), 0.14_dp, 1.0e-12_dp, "cumulative high return")

  grouped(:, 1) = result%bucket_codes(:, 1)
  grouped(:, 2) = 3 - grouped(:, 1)
  returns2(:, 1) = returns(:, 1)
  returns2(:, 2) = 2.0_dp * returns(:, 1)
  config%natural = .false.
  config%by_period = .false.
  call run_grouped_backtest(grouped, returns2, config, grouped_result, status, &
                            by_group=period)
  call assert_true(status == 0, "grouped backtest status")
  call assert_true(grouped_result%n_signals == 2 .and. grouped_result%n_returns == 2, &
                   "multiple signal and return dimensions")
  call assert_close(grouped_result%means(2, 1, 1, 2), 0.09_dp, 1.0e-12_dp, &
                    "multiple-return grouped mean")
  call assert_close(grouped_result%means(1, 2, 1, 1), 0.045_dp, 1.0e-12_dp, &
                    "categorical bucket reversal")

  config%n_buckets = 2
  config%n_by_buckets = 2
  config%by_period = .false.
  config%natural = .false.
  call run_numeric_backtest(signals, returns, config, by_result, status, &
                            by_values=[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                                       5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, &
                                       9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp])
  call assert_true(status == 0 .and. by_result%n_by == 2, "numeric by-variable bucketing")
  call assert_true(sum(by_result%counts(1, 1, 1, :)) == 6 .and. &
                   sum(by_result%counts(1, 1, 2, :)) == 6, "by-variable group counts")

  universe = .true.
  universe(1:4) = .false.
  config%n_by_buckets = 1
  call run_numeric_backtest(signals, returns, config, universe_result, status, &
                            universe=universe)
  call assert_true(status == 0, "universe-filtered backtest status")
  call assert_true(sum(universe_result%counts(1, 1, 1, :)) == 8, "universe filtering")


  call worst_drawdown([0.10_dp, 0.05_dp, 0.08_dp], drawdown, draw_start, draw_end, return_loss)
  call assert_close(drawdown, 1.0_dp - 1.05_dp / 1.10_dp, 1.0e-12_dp, "worst drawdown")
  call assert_true(draw_start == 1 .and. draw_end == 2, "drawdown indices")
  call assert_close(return_loss, -0.05_dp, 1.0e-12_dp, "drawdown return loss")

  print '(a)', "Backtest engine, summaries, spreads, and accessor tests passed."

contains

  pure real(dp) function mean3(x) result(value)
    real(dp), intent(in) :: x(3)
    value = sum(x) / 3.0_dp
  end function mean3

  pure real(dp) function sd3(x) result(value)
    real(dp), intent(in) :: x(3)
    real(dp) :: avg
    avg = mean3(x)
    value = sqrt(sum((x - avg)**2) / 2.0_dp)
  end function sd3

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    call assert_true(abs(actual - expected) <= tolerance, message)
  end subroutine assert_close

  subroutine assert_array_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    call assert_true(size(actual) == size(expected), message // " size")
    call assert_true(all(abs(actual - expected) <= tolerance), message)
  end subroutine assert_array_close

end program test_engine
