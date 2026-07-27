! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

program backtest_csv
  use backtest_kinds, only : dp
  use backtest_engine, only : backtest_config, backtest_result, run_numeric_backtest
  use backtest_summary, only : spread_statistics
  implicit none

  character(len=1024) :: filename, argument, line
  integer :: unit, ios, n, i, status, n_buckets, overlaps
  integer, allocatable :: period(:), ids(:)
  real(dp), allocatable :: signals(:, :), returns(:, :)
  type(backtest_config) :: config
  type(backtest_result) :: result
  real(dp), allocatable :: spread(:), ci_low(:), ci_high(:)

  call get_command_argument(1, filename)
  if (len_trim(filename) == 0) then
    print '(a)', "usage: backtest_csv FILE [N_BUCKETS] [OVERLAPS]"
    print '(a)', "CSV columns: period,id,signal,return"
    stop 2
  end if
  n_buckets = 5
  overlaps = 1
  call get_command_argument(2, argument)
  if (len_trim(argument) > 0) read(argument, *, iostat=ios) n_buckets
  call get_command_argument(3, argument)
  if (len_trim(argument) > 0) read(argument, *, iostat=ios) overlaps

  open(newunit=unit, file=trim(filename), status="old", action="read", iostat=ios)
  if (ios /= 0) error stop "cannot open CSV file"
  read(unit, '(a)', iostat=ios) line
  n = 0
  do
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) > 0) n = n + 1
  end do
  rewind(unit)
  read(unit, '(a)') line
  allocate(period(n), ids(n), signals(n, 1), returns(n, 1))
  do i = 1, n
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) error stop "unexpected end of CSV file"
    call commas_to_spaces(line)
    read(line, *, iostat=ios) period(i), ids(i), signals(i, 1), returns(i, 1)
    if (ios /= 0) error stop "invalid CSV row"
  end do
  close(unit)

  config%n_buckets = n_buckets
  config%by_period = .true.
  config%natural = .true.
  config%overlaps = overlaps
  call run_numeric_backtest(signals, returns, config, result, status, &
                            period=period, ids=ids)
  if (status /= 0) then
    write(*, '(a,i0)') "backtest failed with status ", status
    stop 1
  end if
  call spread_statistics(result, 1, 1, spread, ci_low, ci_high)

  print '(a)', "period,low_mean,high_mean,spread,ci_low,ci_high,turnover"
  do i = 1, result%n_by
    write(*, '(i0,6(",",es14.6))') result%by_levels(i), &
      result%means(1, 1, i, 1), result%means(1, 1, i, n_buckets), &
      spread(i), ci_low(i), ci_high(i), result%turnover(i, 1)
  end do

contains

  subroutine commas_to_spaces(text)
    character(len=*), intent(inout) :: text
    integer :: j
    do j = 1, len_trim(text)
      if (text(j:j) == ',') text(j:j) = ' '
    end do
  end subroutine commas_to_spaces

end program backtest_csv
