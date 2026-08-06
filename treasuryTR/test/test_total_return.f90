! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

program test_total_return
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use treasurytr, only : convexity, dp, mod_duration, period_total_return, total_return, &
    tt_success
  implicit none

  real(dp), parameter :: tol = 2.0e-14_dp
  real(dp) :: yields(4), returns(4), returns_supplied(4)
  real(dp) :: durations(4), convexities(4), expected
  integer :: status
  character(len=160) :: message

  yields = [0.0400_dp, 0.0410_dp, 0.0390_dp, 0.0395_dp]
  durations = mod_duration(yields, 10.0_dp)
  convexities = convexity(yields, 10.0_dp)

  returns = total_return(yields, 10.0_dp, scale=12.0_dp, status=status, message=message)
  call assert_true(status == tt_success, 'total_return status')
  call assert_true(ieee_is_nan(returns(1)), 'first return must be NaN')

  expected = -durations(2) * (yields(2) - yields(1)) + &
    0.5_dp * convexities(2) * (yields(2) - yields(1)) ** 2 + &
    (1.0_dp + yields(1)) ** (1.0_dp / 12.0_dp) - 1.0_dp
  call assert_close(returns(2), expected, tol, 'second return')

  returns_supplied = total_return(yields, 10.0_dp, scale=12.0_dp, mdur=durations, &
    convex=convexities, status=status, message=message)
  call assert_true(status == tt_success, 'supplied sensitivity status')
  call assert_close(maxval(abs(returns(2:) - returns_supplied(2:))), 0.0_dp, tol, &
    'supplied duration and convexity')

  expected = period_total_return(yields(4), yields(3), 10.0_dp, 12.0_dp, &
    source_compatible=.false.)
  returns = total_return(yields, 10.0_dp, scale=12.0_dp, source_compatible=.false.)
  call assert_close(returns(4), expected, tol, 'corrected period return')

  print '(a)', 'test_total_return: PASS'

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

end program test_total_return
