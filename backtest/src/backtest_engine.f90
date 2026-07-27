! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

module backtest_engine
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use backtest_kinds, only : dp
  use backtest_math, only : nan_dp, mean_finite, median_finite, sample_sd_finite, &
                           dense_integer_groups
  use backtest_bucket, only : categorize_quantiles, categorize_by_period, &
                              bucketize_statistics, trim_mask
  use backtest_portfolio, only : turnover_series, overlapping_weights
  implicit none
  private

  type, public :: backtest_config
    integer :: n_buckets = 5
    integer :: n_by_buckets = 1
    logical :: by_period = .true.
    logical :: natural = .false.
    logical :: do_spread = .true.
    integer :: overlaps = 1
    real(dp) :: trim_lower = 0.0025_dp
    real(dp) :: trim_upper = 0.9975_dp
  end type backtest_config

  type, public :: backtest_result
    integer :: n_obs = 0
    integer :: n_signals = 0
    integer :: n_returns = 0
    integer :: n_by = 0
    integer :: n_buckets = 0
    logical :: natural = .false.
    logical :: do_spread = .true.
    logical :: by_period = .false.
    integer :: overlaps = 1
    integer, allocatable :: bucket_codes(:, :)
    integer, allocatable :: by_codes(:)
    integer, allocatable :: by_levels(:)
    real(dp), allocatable :: observation_weights(:)
    real(dp), allocatable :: means(:, :, :, :)
    integer, allocatable :: counts(:, :, :, :)
    real(dp), allocatable :: trimmed_means(:, :, :, :)
    integer, allocatable :: na_counts(:, :, :, :)
    real(dp), allocatable :: return_stats(:, :)
    real(dp), allocatable :: turnover(:, :)
  end type backtest_result

  public :: run_numeric_backtest, run_grouped_backtest

contains

  subroutine run_numeric_backtest(signals, returns, config, result, status, by_group, &
                                  by_values, period, ids, universe)
    real(dp), intent(in) :: signals(:, :), returns(:, :)
    type(backtest_config), intent(in) :: config
    type(backtest_result), intent(out) :: result
    integer, intent(out) :: status
    integer, intent(in), optional :: by_group(:), period(:), ids(:)
    real(dp), intent(in), optional :: by_values(:)
    logical, intent(in), optional :: universe(:)
    real(dp), allocatable :: work_signals(:, :), work_returns(:, :), work_by_values(:), weights(:)
    integer, allocatable :: buckets(:, :), groups(:), local_bucket(:)
    logical, allocatable :: active(:)
    integer :: i, local_status

    status = 0
    if (size(signals, 1) /= size(returns, 1) .or. size(signals, 1) == 0) then
      status = 1
      return
    end if
    if (config%n_buckets < 1 .or. config%n_by_buckets < 1 .or. config%overlaps < 1) then
      status = 2
      return
    end if
    if (merge(1, 0, present(by_group)) + merge(1, 0, present(by_values)) + &
        merge(1, 0, present(period)) > 1) then
      status = 3
      return
    end if
    if (config%by_period .and. .not. present(period)) then
      status = 4
      return
    end if
    if ((config%natural .or. config%overlaps > 1) .and. &
        (.not. present(period) .or. .not. present(ids))) then
      status = 5
      return
    end if
    if (config%overlaps > 1 .and. size(signals, 2) /= 1) then
      status = 6
      return
    end if

    allocate(work_signals, source=signals)
    allocate(work_returns, source=returns)
    allocate(active(size(signals, 1)))
    active = .true.
    if (present(universe)) then
      if (size(universe) /= size(active)) then
        status = 7
        return
      end if
      active = universe
    end if
    do i = 1, size(active)
      if (.not. active(i)) then
        work_signals(i, :) = nan_dp()
        work_returns(i, :) = nan_dp()
      end if
    end do
    if (present(by_values)) then
      allocate(work_by_values, source=by_values)
      do i = 1, size(active)
        if (.not. active(i)) work_by_values(i) = nan_dp()
      end do
    end if

    allocate(buckets(size(signals, 1), size(signals, 2)))
    do i = 1, size(signals, 2)
      if (config%by_period) then
        call categorize_by_period(work_signals(:, i), period, config%n_buckets, &
                                  local_bucket, local_status)
      else
        call categorize_quantiles(work_signals(:, i), config%n_buckets, &
                                  local_bucket, local_status)
      end if
      if (local_status /= 0) then
        status = 10 + local_status
        return
      end if
      buckets(:, i) = local_bucket
      deallocate(local_bucket)
    end do

    if (present(period)) then
      groups = period
    else if (present(by_group)) then
      groups = by_group
    else if (present(by_values)) then
      if (size(by_values) /= size(signals, 1)) then
        status = 8
        return
      end if
      call categorize_quantiles(work_by_values, config%n_by_buckets, local_bucket, local_status)
      if (local_status /= 0) then
        status = 30 + local_status
        return
      end if
      groups = local_bucket
      deallocate(local_bucket)
    else
      allocate(groups(size(signals, 1)))
      groups = 1
    end if

    if (config%overlaps > 1) then
      call overlapping_weights(ids, buckets(:, 1), period, config%n_buckets, &
                               config%overlaps, weights, local_status)
      if (local_status /= 0) then
        status = 20 + local_status
        return
      end if
    else
      allocate(weights(size(signals, 1)))
      weights = 1.0_dp
    end if

    call run_grouped_backtest(buckets, work_returns, config, result, status, &
                              by_group=groups, weights=weights, period=period, ids=ids)
  end subroutine run_numeric_backtest

  subroutine run_grouped_backtest(bucket_codes, returns, config, result, status, &
                                  by_group, weights, period, ids)
    integer, intent(in) :: bucket_codes(:, :)
    real(dp), intent(in) :: returns(:, :)
    type(backtest_config), intent(in) :: config
    type(backtest_result), intent(out) :: result
    integer, intent(out) :: status
    integer, intent(in), optional :: by_group(:), period(:), ids(:)
    real(dp), intent(in), optional :: weights(:)
    integer, allocatable :: groups_raw(:), groups(:), group_levels(:)
    real(dp), allocatable :: weighted_return(:), trim_values(:), local_means(:, :)
    real(dp), allocatable :: local_trim_means(:, :)
    real(dp), allocatable :: obs_weights(:), local_turnover(:)
    integer, allocatable :: local_counts(:, :), local_na(:, :), trim_x(:), trim_y(:)
    logical, allocatable :: keep(:)
    integer :: nobs, nsig, nret, nby, nb, i, j, r, ntrim

    status = 0
    nobs = size(bucket_codes, 1)
    nsig = size(bucket_codes, 2)
    nret = size(returns, 2)
    nb = config%n_buckets
    if (size(returns, 1) /= nobs .or. nobs == 0 .or. nsig == 0 .or. nret == 0) then
      status = 1
      return
    end if
    if (any(bucket_codes < 0) .or. any(bucket_codes > nb)) then
      status = 2
      return
    end if
    allocate(groups_raw(nobs))
    if (present(by_group)) then
      if (size(by_group) /= nobs) then
        status = 3
        return
      end if
      groups_raw = by_group
    else
      groups_raw = 1
    end if
    call dense_integer_groups(groups_raw, groups, group_levels)
    nby = size(group_levels)

    allocate(obs_weights(nobs))
    if (present(weights)) then
      if (size(weights) /= nobs) then
        status = 4
        return
      end if
      obs_weights = weights
    else
      obs_weights = 1.0_dp
    end if

    result%n_obs = nobs
    result%n_signals = nsig
    result%n_returns = nret
    result%n_by = nby
    result%n_buckets = nb
    result%natural = config%natural
    result%do_spread = config%do_spread
    result%by_period = config%by_period
    result%overlaps = config%overlaps
    allocate(result%bucket_codes, source=bucket_codes)
    allocate(result%by_codes, source=groups)
    allocate(result%by_levels, source=group_levels)
    allocate(result%observation_weights, source=obs_weights)
    allocate(result%means(nret, nsig, nby, nb))
    allocate(result%counts(nret, nsig, nby, nb))
    allocate(result%trimmed_means(nret, nsig, nby, nb))
    allocate(result%na_counts(nret, nsig, nby, nb))
    allocate(result%return_stats(nret, 6))
    result%means = nan_dp()
    result%counts = 0
    result%trimmed_means = nan_dp()
    result%na_counts = 0
    result%return_stats = nan_dp()

    allocate(weighted_return(nobs))
    do r = 1, nret
      do i = 1, nobs
        if (ieee_is_finite(returns(i, r)) .and. ieee_is_finite(obs_weights(i))) then
          weighted_return(i) = returns(i, r) * obs_weights(i)
        else
          weighted_return(i) = nan_dp()
        end if
      end do
      if (count(ieee_is_finite(weighted_return)) > 0) then
        result%return_stats(r, 1) = minval(weighted_return, mask=ieee_is_finite(weighted_return))
        result%return_stats(r, 2) = maxval(weighted_return, mask=ieee_is_finite(weighted_return))
      else
        result%return_stats(r, 1:2) = nan_dp()
      end if
      result%return_stats(r, 3) = mean_finite(weighted_return)
      result%return_stats(r, 4) = median_finite(weighted_return)
      result%return_stats(r, 5) = sample_sd_finite(weighted_return)
      result%return_stats(r, 6) = real(count(.not. ieee_is_finite(weighted_return)), dp)
      call trim_mask(weighted_return, config%trim_lower, config%trim_upper, keep)
      ntrim = count(keep)
      allocate(trim_values(ntrim), trim_x(ntrim), trim_y(ntrim))
      do j = 1, nsig
        allocate(local_means(nby, nb), local_counts(nby, nb), local_na(nby, nb))
        call bucketize_statistics(weighted_return, bucket_codes(:, j), groups, nb, nby, &
                                  local_means, local_counts, local_na)
        result%means(r, j, :, :) = local_means
        result%counts(r, j, :, :) = local_counts
        result%na_counts(r, j, :, :) = local_na
        i = 0
        do ntrim = 1, nobs
          if (keep(ntrim)) then
            i = i + 1
            trim_values(i) = weighted_return(ntrim)
            trim_x(i) = bucket_codes(ntrim, j)
            trim_y(i) = groups(ntrim)
          end if
        end do
        allocate(local_trim_means(nby, nb))
        local_counts = 0
        local_na = 0
        call bucketize_statistics(trim_values, trim_x, trim_y, nb, nby, &
                                  local_trim_means, local_counts, local_na)
        result%trimmed_means(r, j, :, :) = local_trim_means
        deallocate(local_means, local_counts, local_na, local_trim_means)
      end do
      deallocate(keep, trim_values, trim_x, trim_y)
    end do

    if (config%natural) then
      if (.not. present(period) .or. .not. present(ids)) then
        status = 5
        return
      end if
      allocate(result%turnover(nby, nsig))
      result%turnover = nan_dp()
      do j = 1, nsig
        call turnover_series(ids, bucket_codes(:, j), period, nb, local_turnover)
        if (size(local_turnover) /= nby) then
          status = 6
          return
        end if
        result%turnover(:, j) = local_turnover
        deallocate(local_turnover)
      end do
    else
      allocate(result%turnover(0, 0))
    end if
  end subroutine run_grouped_backtest

end module backtest_engine
