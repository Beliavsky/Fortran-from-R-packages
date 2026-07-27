! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

program test_portfolio
  use backtest_kinds, only : dp
  use backtest_portfolio, only : turnover_series, tribucket_weights, &
                                 scale_period_weights, overlapping_weights
  implicit none

  integer :: ids(12), buckets(12), periods(12), status
  real(dp), allocatable :: turnover(:), base(:), scaled(:), overlap_weight(:)
  real(dp) :: raw(8)

  ids = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]
  periods = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
  buckets = [1, 1, 3, 3, 1, 3, 1, 3, 3, 1, 1, 3]

  call turnover_series(ids, buckets, periods, 3, turnover)
  call assert_true(size(turnover) == 3, "turnover periods")
  call assert_close(turnover(2), 0.5_dp, 1.0e-12_dp, "turnover period 2")
  call assert_close(turnover(3), 0.5_dp, 1.0e-12_dp, "turnover period 3")

  call tribucket_weights(buckets, periods, 3, base)
  call assert_close(base(1), -0.5_dp, 1.0e-12_dp, "low base weight")
  call assert_close(base(3), 0.5_dp, 1.0e-12_dp, "high base weight")
  call assert_close(base(6), 0.5_dp, 1.0e-12_dp, "second-period high weight")

  raw = [-2.0_dp, -1.0_dp, 1.0_dp, 3.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp]
  call scale_period_weights(raw, [1, 1, 1, 1, 2, 2, 2, 2], scaled)
  call assert_close(scaled(1), 2.0_dp / 3.0_dp, 1.0e-12_dp, "negative scaling 1")
  call assert_close(scaled(2), 1.0_dp / 3.0_dp, 1.0e-12_dp, "negative scaling 2")
  call assert_close(scaled(3), 0.25_dp, 1.0e-12_dp, "positive scaling 1")
  call assert_close(scaled(4), 0.75_dp, 1.0e-12_dp, "positive scaling 2")
  call assert_close(scaled(5), 1.0_dp, 1.0e-12_dp, "single negative scaling")
  call assert_close(scaled(6), 0.5_dp, 1.0e-12_dp, "zero scaling 1")
  call assert_close(scaled(7), 0.5_dp, 1.0e-12_dp, "zero scaling 2")
  call assert_close(scaled(8), 1.0_dp, 1.0e-12_dp, "single positive scaling")

  call overlapping_weights(ids, buckets, periods, 3, 2, overlap_weight, status)
  call assert_true(status == 0, "overlap weight status")
  call assert_array_close(overlap_weight(1:4), [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], &
                          1.0e-12_dp, "first-period overlap weights")
  call assert_array_close(overlap_weight(5:8), [1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], &
                          1.0e-12_dp, "second-period overlap weights")
  call assert_array_close(overlap_weight(9:12), [0.5_dp, 0.5_dp, 1.0_dp, 1.0_dp], &
                          1.0e-12_dp, "third-period overlap weights")

  print '(a)', "Turnover and overlapping-portfolio tests passed."

contains

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

end program test_portfolio
