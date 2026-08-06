! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

program test_preprocessing
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
  use treasurytr, only : carry_forward, dp, percent_to_decimal, prepare_yields
  implicit none

  real(dp), parameter :: tol = 1.0e-15_dp
  real(dp) :: raw(5), prepared(5), expected(5), carried(5), nan_value
  real(dp) :: matrix_raw(4, 2), matrix_prepared(4, 2)

  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  raw = [nan_value, 4.0_dp, nan_value, nan_value, 4.2_dp]
  prepared = prepare_yields(raw)
  expected = [nan_value, 0.04_dp, 0.04_dp, 0.04_dp, 0.042_dp]

  call assert_true(ieee_is_nan(prepared(1)), 'leading missing value')
  call assert_close(maxval(abs(prepared(2:) - expected(2:))), 0.0_dp, tol, &
    'LOCF and percent adjustment')
  call assert_close(percent_to_decimal(5.0_dp), 0.05_dp, tol, 'percent conversion')
  carried = carry_forward(raw)
  call assert_close(carried(5), 4.2_dp, tol, 'vector carry forward')

  matrix_raw(:, 1) = [3.0_dp, nan_value, 3.2_dp, nan_value]
  matrix_raw(:, 2) = [4.0_dp, 4.1_dp, nan_value, 4.3_dp]
  matrix_prepared = prepare_yields(matrix_raw)
  call assert_close(matrix_prepared(2, 1), 0.03_dp, tol, 'matrix LOCF column one')
  call assert_close(matrix_prepared(3, 2), 0.041_dp, tol, 'matrix LOCF column two')

  print '(a)', 'test_preprocessing: PASS'

contains

  subroutine assert_close(actual, expected_value, tolerance, label)
    real(dp), intent(in) :: actual, expected_value, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected_value) > tolerance) then
      print '(a,2(1x,es24.16))', trim(label), actual, expected_value
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      print '(a)', trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_preprocessing
