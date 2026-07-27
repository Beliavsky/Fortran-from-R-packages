! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

program demo_backtest
  use backtest_kinds, only : dp
  use backtest_engine, only : backtest_config, backtest_result, run_numeric_backtest
  use backtest_summary, only : spread_statistics, sharpe_ratios
  implicit none

  real(dp) :: signal(20, 1), returns(20, 1)
  integer :: period(20), ids(20), i, p, status
  type(backtest_config) :: config
  type(backtest_result) :: result
  real(dp), allocatable :: spread(:), ci_low(:), ci_high(:), sharpe(:)

  do p = 1, 5
    do i = 1, 4
      period((p - 1) * 4 + i) = 202600 + p
      ids((p - 1) * 4 + i) = i
      signal((p - 1) * 4 + i, 1) = real(mod(i + p, 4) + 1, dp)
      returns((p - 1) * 4 + i, 1) = 0.002_dp * real(i, dp) + &
                                             0.001_dp * real(p, dp)
    end do
  end do

  config%n_buckets = 2
  config%by_period = .true.
  config%natural = .true.
  call run_numeric_backtest(signal, returns, config, result, status, &
                            period=period, ids=ids)
  if (status /= 0) error stop "demo backtest failed"

  call spread_statistics(result, 1, 1, spread, ci_low, ci_high)
  call sharpe_ratios(result, 1, sharpe)

  print '(a)', "period     low mean     high mean    spread       turnover"
  do p = 1, result%n_by
    write(*, '(i8,4(1x,f12.6))') result%by_levels(p), &
      result%means(1, 1, p, 1), result%means(1, 1, p, 2), spread(p), &
      result%turnover(p, 1)
  end do
  write(*, '(a,f12.6)') "raw Sharpe ratio: ", sharpe(1)

end program demo_backtest
