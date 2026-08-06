! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

program test_duration_convexity
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use treasurytr, only : convexity, dp, mod_duration
  implicit none

  real(dp), parameter :: tol = 1.0e-13_dp
  real(dp) :: y, m, expected_duration, expected_convexity, z
  real(dp) :: values(3)

  y = 0.04_dp
  m = 10.0_dp
  z = 1.0_dp + 0.5_dp * y
  expected_duration = (1.0_dp - z ** (-2.0_dp * m)) / y
  expected_convexity = 2.0_dp / y ** 2 * (1.0_dp - z ** (-2.0_dp * m)) - &
    2.0_dp * m / y * z ** (-2.0_dp * m - 1.0_dp)

  call assert_close(mod_duration(y, m), expected_duration, tol, 'duration formula')
  call assert_close(convexity(y, m), expected_convexity, tol, 'convexity formula')
  call assert_close(mod_duration(0.0_dp, m, .false.), m, tol, 'zero-yield duration limit')
  call assert_close(convexity(0.0_dp, m, .false.), 105.0_dp, tol, &
    'zero-yield convexity limit')
  call assert_true(ieee_is_nan(mod_duration(0.0_dp, m)), &
    'source-compatible zero duration must be NaN')
  call assert_true(ieee_is_nan(convexity(0.0_dp, m)), &
    'source-compatible zero convexity must be NaN')

  values = mod_duration([0.02_dp, 0.04_dp, 0.06_dp], 5.0_dp)
  call assert_true(all(values > 0.0_dp), 'elemental duration evaluation')
  call assert_close(mod_duration(1.0e-10_dp, m, .false.), m, 1.0e-8_dp, &
    'stable near-zero duration')
  call assert_close(convexity(1.0e-10_dp, m, .false.), 105.0_dp, 1.0e-7_dp, &
    'stable near-zero convexity')

  print '(a)', 'test_duration_convexity: PASS'

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      print '(a,2(1x,es24.16))', trim(label), actual, expected
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

end program test_duration_convexity
