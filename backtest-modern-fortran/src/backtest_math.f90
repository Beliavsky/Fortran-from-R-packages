! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

module backtest_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use backtest_kinds, only : dp
  implicit none
  private

  public :: nan_dp, mean_finite, sample_sd_finite, median_finite
  public :: quantile_type7, sort_real, dense_integer_groups
  public :: cumulative_sum_columns, worst_drawdown

contains

  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key

    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_real

  pure real(dp) function mean_finite(x) result(value)
    real(dp), intent(in) :: x(:)
    integer :: i, n
    real(dp) :: total

    total = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        total = total + x(i)
        n = n + 1
      end if
    end do
    if (n > 0) then
      value = total / real(n, dp)
    else
      value = nan_dp()
    end if
  end function mean_finite

  pure real(dp) function sample_sd_finite(x) result(value)
    real(dp), intent(in) :: x(:)
    integer :: i, n
    real(dp) :: avg, ss

    avg = mean_finite(x)
    if (.not. ieee_is_finite(avg)) then
      value = nan_dp()
      return
    end if
    n = 0
    ss = 0.0_dp
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        ss = ss + (x(i) - avg)**2
        n = n + 1
      end if
    end do
    if (n > 1) then
      value = sqrt(ss / real(n - 1, dp))
    else
      value = nan_dp()
    end if
  end function sample_sd_finite

  real(dp) function quantile_type7(x, p) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: p
    real(dp), allocatable :: y(:)
    integer :: i, n, j
    real(dp) :: h, g

    n = count(ieee_is_finite(x))
    if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = nan_dp()
      return
    end if
    allocate(y(n))
    j = 0
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        j = j + 1
        y(j) = x(i)
      end if
    end do
    call sort_real(y)
    if (n == 1) then
      value = y(1)
      return
    end if
    h = 1.0_dp + real(n - 1, dp) * p
    j = int(floor(h))
    g = h - real(j, dp)
    if (j >= n) then
      value = y(n)
    else
      value = (1.0_dp - g) * y(j) + g * y(j + 1)
    end if
  end function quantile_type7

  real(dp) function median_finite(x) result(value)
    real(dp), intent(in) :: x(:)
    value = quantile_type7(x, 0.5_dp)
  end function median_finite

  subroutine dense_integer_groups(raw, dense, levels)
    integer, intent(in) :: raw(:)
    integer, allocatable, intent(out) :: dense(:)
    integer, allocatable, intent(out) :: levels(:)
    integer, allocatable :: work(:)
    integer :: i, j, nlev, key

    allocate(work(size(raw)))
    work = raw
    do i = 2, size(work)
      key = work(i)
      j = i - 1
      do while (j >= 1)
        if (work(j) <= key) exit
        work(j + 1) = work(j)
        j = j - 1
      end do
      work(j + 1) = key
    end do
    nlev = 0
    do i = 1, size(work)
      if (i == 1) then
        nlev = nlev + 1
      else if (work(i) /= work(i - 1)) then
        nlev = nlev + 1
      end if
    end do
    allocate(levels(nlev), dense(size(raw)))
    j = 0
    do i = 1, size(work)
      if (i == 1) then
        j = j + 1
        levels(j) = work(i)
      else if (work(i) /= work(i - 1)) then
        j = j + 1
        levels(j) = work(i)
      end if
    end do
    do i = 1, size(raw)
      do j = 1, nlev
        if (raw(i) == levels(j)) then
          dense(i) = j
          exit
        end if
      end do
    end do
  end subroutine dense_integer_groups

  subroutine cumulative_sum_columns(x, cumulative)
    real(dp), intent(in) :: x(:, :)
    real(dp), allocatable, intent(out) :: cumulative(:, :)
    integer :: i, j

    allocate(cumulative(size(x, 1), size(x, 2)))
    cumulative = nan_dp()
    do j = 1, size(x, 2)
      if (size(x, 1) > 0) cumulative(1, j) = x(1, j)
      do i = 2, size(x, 1)
        if (ieee_is_finite(x(i, j)) .and. ieee_is_finite(cumulative(i - 1, j))) then
          cumulative(i, j) = cumulative(i - 1, j) + x(i, j)
        else if (ieee_is_finite(x(i, j))) then
          cumulative(i, j) = x(i, j)
        else
          cumulative(i, j) = cumulative(i - 1, j)
        end if
      end do
    end do
  end subroutine cumulative_sum_columns

  subroutine worst_drawdown(cumulative_return, drawdown, start_index, end_index, return_loss)
    real(dp), intent(in) :: cumulative_return(:)
    real(dp), intent(out) :: drawdown
    integer, intent(out) :: start_index, end_index
    real(dp), intent(out) :: return_loss
    real(dp) :: wealth, peak, current_dd
    integer :: i, peak_index

    drawdown = 0.0_dp
    return_loss = 0.0_dp
    start_index = 1
    end_index = 1
    peak = 1.0_dp
    peak_index = 1
    do i = 1, size(cumulative_return)
      if (.not. ieee_is_finite(cumulative_return(i))) cycle
      wealth = 1.0_dp + cumulative_return(i)
      if (wealth >= peak) then
        peak = wealth
        peak_index = i
      end if
      if (peak > 0.0_dp) then
        current_dd = 1.0_dp - wealth / peak
        if (current_dd >= drawdown) then
          drawdown = current_dd
          start_index = peak_index
          end_index = i
          return_loss = cumulative_return(i) - cumulative_return(peak_index)
        end if
      end if
    end do
  end subroutine worst_drawdown

end module backtest_math
