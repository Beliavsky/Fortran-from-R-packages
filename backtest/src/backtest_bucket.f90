! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

module backtest_bucket
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use backtest_kinds, only : dp
  use backtest_math, only : nan_dp, quantile_type7, dense_integer_groups
  implicit none
  private

  public :: categorize_quantiles, categorize_by_period
  public :: bucketize_statistics, trim_mask

contains

  subroutine categorize_quantiles(x, n_buckets, bucket, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n_buckets
    integer, allocatable, intent(out) :: bucket(:)
    integer, intent(out) :: status
    real(dp), allocatable :: breaks(:)
    integer :: i, j

    status = 0
    allocate(bucket(size(x)))
    bucket = 0
    if (n_buckets < 1 .or. count(ieee_is_finite(x)) < n_buckets) then
      status = 1
      return
    end if
    allocate(breaks(0:n_buckets))
    do j = 0, n_buckets
      breaks(j) = quantile_type7(x, real(j, dp) / real(n_buckets, dp))
    end do
    do j = 1, n_buckets
      if (breaks(j) <= breaks(j - 1)) then
        status = 2
        return
      end if
    end do
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) cycle
      do j = 1, n_buckets
        if (x(i) <= breaks(j)) then
          bucket(i) = j
          exit
        end if
      end do
      if (bucket(i) == 0) bucket(i) = n_buckets
    end do
  end subroutine categorize_quantiles

  subroutine categorize_by_period(x, period, n_buckets, bucket, status)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: period(:)
    integer, intent(in) :: n_buckets
    integer, allocatable, intent(out) :: bucket(:)
    integer, intent(out) :: status
    integer, allocatable :: dense(:), levels(:), index(:), local_bucket(:)
    real(dp), allocatable :: local_x(:)
    integer :: g, i, n, local_status

    status = 0
    if (size(x) /= size(period)) then
      status = 3
      allocate(bucket(0))
      return
    end if
    call dense_integer_groups(period, dense, levels)
    allocate(bucket(size(x)))
    bucket = 0
    do g = 1, size(levels)
      n = count(dense == g)
      allocate(local_x(n), index(n))
      n = 0
      do i = 1, size(x)
        if (dense(i) == g) then
          n = n + 1
          local_x(n) = x(i)
          index(n) = i
        end if
      end do
      call categorize_quantiles(local_x, n_buckets, local_bucket, local_status)
      if (local_status /= 0) then
        status = local_status
        return
      end if
      do i = 1, size(index)
        bucket(index(i)) = local_bucket(i)
      end do
      deallocate(local_x, index, local_bucket)
    end do
  end subroutine categorize_by_period

  subroutine bucketize_statistics(values, x_group, y_group, n_x, n_y, means, counts, na_counts)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: x_group(:), y_group(:)
    integer, intent(in) :: n_x, n_y
    real(dp), intent(out) :: means(n_y, n_x)
    integer, intent(out) :: counts(n_y, n_x), na_counts(n_y, n_x)
    real(dp) :: totals(n_y, n_x)
    integer :: finite_counts(n_y, n_x)
    integer :: i, xb, yb

    totals = 0.0_dp
    finite_counts = 0
    counts = 0
    na_counts = 0
    means = nan_dp()
    do i = 1, size(values)
      xb = x_group(i)
      yb = y_group(i)
      if (xb < 1 .or. xb > n_x .or. yb < 1 .or. yb > n_y) cycle
      counts(yb, xb) = counts(yb, xb) + 1
      if (ieee_is_finite(values(i))) then
        totals(yb, xb) = totals(yb, xb) + values(i)
        finite_counts(yb, xb) = finite_counts(yb, xb) + 1
      else
        na_counts(yb, xb) = na_counts(yb, xb) + 1
      end if
    end do
    do yb = 1, n_y
      do xb = 1, n_x
        if (finite_counts(yb, xb) > 0) then
          means(yb, xb) = totals(yb, xb) / real(finite_counts(yb, xb), dp)
        end if
      end do
    end do
  end subroutine bucketize_statistics

  subroutine trim_mask(x, lower_probability, upper_probability, keep)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: lower_probability, upper_probability
    logical, allocatable, intent(out) :: keep(:)
    real(dp) :: lower, upper
    integer :: i

    lower = quantile_type7(x, lower_probability)
    upper = quantile_type7(x, upper_probability)
    allocate(keep(size(x)))
    do i = 1, size(x)
      keep(i) = ieee_is_finite(x(i)) .and. x(i) > lower .and. x(i) < upper
    end do
  end subroutine trim_mask

end module backtest_bucket
