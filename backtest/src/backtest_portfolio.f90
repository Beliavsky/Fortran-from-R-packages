! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

module backtest_portfolio
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use backtest_kinds, only : dp
  use backtest_math, only : nan_dp, dense_integer_groups
  implicit none
  private

  public :: turnover_series, tribucket_weights, scale_period_weights, overlapping_weights

contains

  subroutine turnover_series(ids, buckets, periods, n_buckets, turnover, period_levels)
    integer, intent(in) :: ids(:), buckets(:), periods(:)
    integer, intent(in) :: n_buckets
    real(dp), allocatable, intent(out) :: turnover(:)
    integer, allocatable, intent(out), optional :: period_levels(:)
    integer, allocatable :: dense(:), levels(:)
    integer :: g, side, i, j, id_value, n_old, n_new
    integer, allocatable :: unique_ids(:)
    real(dp) :: old_weight, new_weight, total

    call dense_integer_groups(periods, dense, levels)
    allocate(turnover(size(levels)))
    turnover = nan_dp()
    if (present(period_levels)) then
      allocate(period_levels(size(levels)))
      period_levels = levels
    end if
    do g = 2, size(levels)
      allocate(unique_ids(size(ids) * 2))
      unique_ids = 0
      j = 0
      do side = 1, 2
        do i = 1, size(ids)
          if ((dense(i) == g .or. dense(i) == g - 1) .and. &
              ((side == 1 .and. buckets(i) == n_buckets) .or. &
               (side == 2 .and. buckets(i) == 1))) then
            if (.not. any(unique_ids(1:j) == ids(i))) then
              j = j + 1
              unique_ids(j) = ids(i)
            end if
          end if
        end do
      end do
      total = 0.0_dp
      do side = 1, 2
        n_old = count(dense == g - 1 .and. merge(buckets == n_buckets, buckets == 1, side == 1))
        n_new = count(dense == g .and. merge(buckets == n_buckets, buckets == 1, side == 1))
        do i = 1, j
          id_value = unique_ids(i)
          old_weight = 0.0_dp
          new_weight = 0.0_dp
          if (n_old > 0) then
            if (any(dense == g - 1 .and. ids == id_value .and. &
                    merge(buckets == n_buckets, buckets == 1, side == 1))) then
              old_weight = 1.0_dp / real(n_old, dp)
            end if
          end if
          if (n_new > 0) then
            if (any(dense == g .and. ids == id_value .and. &
                    merge(buckets == n_buckets, buckets == 1, side == 1))) then
              new_weight = 1.0_dp / real(n_new, dp)
            end if
          end if
          total = total + abs(old_weight - new_weight)
        end do
      end do
      turnover(g) = total / 4.0_dp
      deallocate(unique_ids)
    end do
  end subroutine turnover_series

  subroutine tribucket_weights(buckets, periods, n_buckets, weights)
    integer, intent(in) :: buckets(:), periods(:), n_buckets
    real(dp), allocatable, intent(out) :: weights(:)
    integer, allocatable :: dense(:), levels(:)
    integer :: g, i, n_low, n_high

    call dense_integer_groups(periods, dense, levels)
    allocate(weights(size(buckets)))
    weights = 0.0_dp
    do g = 1, size(levels)
      n_low = count(dense == g .and. buckets == 1)
      n_high = count(dense == g .and. buckets == n_buckets)
      do i = 1, size(buckets)
        if (dense(i) /= g) cycle
        if (buckets(i) == 1 .and. n_low > 0) weights(i) = -1.0_dp / real(n_low, dp)
        if (buckets(i) == n_buckets .and. n_high > 0) weights(i) = 1.0_dp / real(n_high, dp)
      end do
    end do
  end subroutine tribucket_weights

  subroutine scale_period_weights(raw_weights, periods, weights)
    real(dp), intent(in) :: raw_weights(:)
    integer, intent(in) :: periods(:)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, allocatable :: dense(:), levels(:)
    integer :: g, i, n_zero
    real(dp) :: positive_sum, negative_sum, tol

    call dense_integer_groups(periods, dense, levels)
    allocate(weights, source=raw_weights)
    tol = epsilon(1.0_dp)
    do g = 1, size(levels)
      positive_sum = sum(weights, mask=dense == g .and. weights > tol)
      negative_sum = sum(weights, mask=dense == g .and. weights < -tol)
      n_zero = count(dense == g .and. abs(weights) <= tol)
      do i = 1, size(weights)
        if (dense(i) /= g) cycle
        if (weights(i) > tol .and. abs(positive_sum) > tol) then
          weights(i) = weights(i) / positive_sum
        else if (weights(i) < -tol .and. abs(negative_sum) > tol) then
          weights(i) = weights(i) / negative_sum
        else if (abs(weights(i)) <= tol .and. n_zero > 0) then
          weights(i) = 1.0_dp / real(n_zero, dp)
        end if
      end do
    end do
  end subroutine scale_period_weights

  subroutine overlapping_weights(ids, buckets, periods, n_buckets, overlaps, weights, status)
    integer, intent(in) :: ids(:), buckets(:), periods(:), n_buckets, overlaps
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    integer, allocatable :: dense(:), levels(:)
    real(dp), allocatable :: base(:), raw(:)
    integer :: i, j, first_group, this_group

    status = 0
    if (size(ids) /= size(buckets) .or. size(periods) /= size(buckets) .or. overlaps < 1) then
      status = 1
      allocate(weights(0))
      return
    end if
    call dense_integer_groups(periods, dense, levels)
    if (overlaps > size(levels)) then
      status = 2
      allocate(weights(0))
      return
    end if
    call tribucket_weights(buckets, periods, n_buckets, base)
    allocate(raw(size(base)))
    raw = 0.0_dp
    do i = 1, size(base)
      this_group = dense(i)
      first_group = max(1, this_group - overlaps + 1)
      do j = 1, size(base)
        if (ids(j) == ids(i) .and. dense(j) >= first_group .and. dense(j) <= this_group) then
          raw(i) = raw(i) + base(j)
        end if
      end do
    end do
    call scale_period_weights(raw, periods, weights)
  end subroutine overlapping_weights

end module backtest_portfolio
