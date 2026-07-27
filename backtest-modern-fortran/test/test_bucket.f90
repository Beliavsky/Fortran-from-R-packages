! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of backtest-modern-fortran.
! It is free software: you may redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 2 of the License, or any later version.

program test_bucket
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use backtest_kinds, only : dp
  use backtest_math, only : nan_dp, quantile_type7
  use backtest_bucket, only : categorize_quantiles, categorize_by_period, &
                              bucketize_statistics, trim_mask
  implicit none

  real(dp) :: x(10), values(6), means(2, 2)
  integer :: periods(8), groups_x(6), groups_y(6)
  integer :: counts(2, 2), na_counts(2, 2), status, i
  integer, allocatable :: bucket(:)
  logical, allocatable :: keep(:)

  x = [(real(i, dp), i = 1, 10)]
  call categorize_quantiles(x, 5, bucket, status)
  call assert_true(status == 0, "quantile categorization status")
  call assert_true(all(bucket == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]), &
                   "quantile bucket assignments")
  call assert_close(quantile_type7(x, 0.2_dp), 2.8_dp, 1.0e-12_dp, &
                    "R type-7 quantile")

  call categorize_quantiles([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 2, bucket, status)
  call assert_true(status == 2, "duplicate quantile breaks")

  periods = [1, 1, 1, 1, 2, 2, 2, 2]
  call categorize_by_period([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                             10.0_dp, 11.0_dp, 12.0_dp, 13.0_dp], &
                            periods, 2, bucket, status)
  call assert_true(status == 0, "period categorization status")
  call assert_true(all(bucket == [1, 1, 2, 2, 1, 1, 2, 2]), &
                   "period bucket assignments")

  values = [1.0_dp, nan_dp(), 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  groups_x = [1, 1, 2, 2, 1, 2]
  groups_y = [1, 1, 1, 1, 2, 2]
  call bucketize_statistics(values, groups_x, groups_y, 2, 2, means, counts, na_counts)
  call assert_close(means(1, 1), 1.0_dp, 1.0e-12_dp, "bucket mean 1")
  call assert_close(means(1, 2), 3.5_dp, 1.0e-12_dp, "bucket mean 2")
  call assert_close(means(2, 1), 5.0_dp, 1.0e-12_dp, "bucket mean 3")
  call assert_close(means(2, 2), 6.0_dp, 1.0e-12_dp, "bucket mean 4")
  call assert_true(counts(1, 1) == 2 .and. na_counts(1, 1) == 1, &
                   "counts include missing values")
  call assert_true(.not. ieee_is_nan(means(1, 1)), "finite grouped mean")

  call trim_mask([(real(i, dp), i = 1, 100)], 0.1_dp, 0.9_dp, keep)
  call assert_true(count(keep) == 80, "strict trimmed range")
  call assert_true(keep(11) .and. keep(90) .and. .not. keep(10) .and. .not. keep(91), &
                   "trim endpoints")

  print '(a)', "Bucket, quantile, grouping, and trimming tests passed."

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

end program test_bucket
