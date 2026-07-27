! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

program overlap_example
  use backtest_kinds, only : dp
  use backtest_engine, only : backtest_config, backtest_result, run_numeric_backtest
  use backtest_summary, only : spread_statistics
  implicit none

  integer :: ids(12), period(12), status, i
  real(dp) :: signal(12, 1), returns(12, 1)
  type(backtest_config) :: config
  type(backtest_result) :: result
  real(dp), allocatable :: spread(:), low(:), high(:)

  ids = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]
  period = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
  signal(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                  1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp, &
                  4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
  do i = 1, 12
    returns(i, 1) = 0.001_dp * real(i, dp)
  end do

  config%n_buckets = 3
  config%by_period = .true.
  config%overlaps = 2
  call run_numeric_backtest(signal, returns, config, result, status, &
                            period=period, ids=ids)
  if (status /= 0) error stop "overlap example failed"
  call spread_statistics(result, 1, 1, spread, low, high)

  print '(a)', "observation weights:"
  write(*, '(*(f8.4,1x))') result%observation_weights
  print '(a)', "period spreads:"
  write(*, '(*(f10.6,1x))') spread

end program overlap_example
