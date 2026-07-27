! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

module backtest_summary
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use backtest_kinds, only : dp
  use backtest_math, only : nan_dp, mean_finite, sample_sd_finite, cumulative_sum_columns
  use backtest_engine, only : backtest_result
  implicit none
  private

  public :: spread_statistics, sharpe_ratios, total_counts
  public :: marginal_counts, mean_rows, cumulative_bucket_returns

contains

  subroutine spread_statistics(result, return_index, signal_index, spread, ci_low, ci_high)
    type(backtest_result), intent(in) :: result
    integer, intent(in) :: return_index, signal_index
    real(dp), allocatable, intent(out) :: spread(:), ci_low(:), ci_high(:)
    integer :: g, total_n
    real(dp) :: se, sigma

    allocate(spread(result%n_by), ci_low(result%n_by), ci_high(result%n_by))
    sigma = result%return_stats(return_index, 5)
    do g = 1, result%n_by
      spread(g) = result%means(return_index, signal_index, g, result%n_buckets) - &
                  result%means(return_index, signal_index, g, 1)
      total_n = result%counts(return_index, signal_index, g, 1) + &
                result%counts(return_index, signal_index, g, result%n_buckets)
      if (total_n > 0 .and. ieee_is_finite(sigma)) then
        se = sigma / sqrt(real(total_n, dp))
        ci_low(g) = spread(g) - 2.0_dp * se
        ci_high(g) = spread(g) + 2.0_dp * se
      else
        ci_low(g) = nan_dp()
        ci_high(g) = nan_dp()
      end if
    end do
  end subroutine spread_statistics

  subroutine sharpe_ratios(result, return_index, sharpe)
    type(backtest_result), intent(in) :: result
    integer, intent(in) :: return_index
    real(dp), allocatable, intent(out) :: sharpe(:)
    real(dp), allocatable :: spreads(:)
    real(dp) :: sd
    integer :: j, g

    allocate(sharpe(result%n_signals), spreads(result%n_by))
    do j = 1, result%n_signals
      do g = 1, result%n_by
        spreads(g) = result%means(return_index, j, g, result%n_buckets) - &
                     result%means(return_index, j, g, 1)
      end do
      sd = sample_sd_finite(spreads)
      if (ieee_is_finite(sd) .and. sd > 0.0_dp) then
        sharpe(j) = mean_finite(spreads) / sd
      else
        sharpe(j) = nan_dp()
      end if
    end do
  end subroutine sharpe_ratios

  subroutine total_counts(result, return_index, low_high_only, totals)
    type(backtest_result), intent(in) :: result
    integer, intent(in) :: return_index
    logical, intent(in) :: low_high_only
    integer, allocatable, intent(out) :: totals(:, :)
    integer :: g, j

    allocate(totals(result%n_by, result%n_signals))
    do j = 1, result%n_signals
      do g = 1, result%n_by
        if (low_high_only) then
          totals(g, j) = result%counts(return_index, j, g, 1) + &
                         result%counts(return_index, j, g, result%n_buckets)
        else
          totals(g, j) = sum(result%counts(return_index, j, g, :))
        end if
      end do
    end do
  end subroutine total_counts

  subroutine marginal_counts(result, return_index, signal_index, table)
    type(backtest_result), intent(in) :: result
    integer, intent(in) :: return_index, signal_index
    integer, allocatable, intent(out) :: table(:, :)
    integer :: g, b

    allocate(table(result%n_by + 1, result%n_buckets + 1))
    table = 0
    table(1:result%n_by, 1:result%n_buckets) = &
      result%counts(return_index, signal_index, :, :)
    do g = 1, result%n_by
      table(g, result%n_buckets + 1) = sum(table(g, 1:result%n_buckets))
    end do
    do b = 1, result%n_buckets + 1
      table(result%n_by + 1, b) = sum(table(1:result%n_by, b))
    end do
  end subroutine marginal_counts

  subroutine mean_rows(x, row_mean)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: row_mean(:)
    integer :: j

    allocate(row_mean(size(x, 2)))
    do j = 1, size(x, 2)
      row_mean(j) = mean_finite(x(:, j))
    end do
  end subroutine mean_rows

  subroutine cumulative_bucket_returns(result, return_index, signal_index, cumulative)
    type(backtest_result), intent(in) :: result
    integer, intent(in) :: return_index, signal_index
    real(dp), allocatable, intent(out) :: cumulative(:, :)
    real(dp), allocatable :: period_means(:, :)

    allocate(period_means(result%n_by, result%n_buckets))
    period_means = result%means(return_index, signal_index, :, :)
    call cumulative_sum_columns(period_means, cumulative)
  end subroutine cumulative_bucket_returns

end module backtest_summary
