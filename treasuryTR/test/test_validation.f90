! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

program test_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use treasurytr, only : convexity, dp, mod_duration, total_return, tt_err_maturity, &
    tt_err_scale, tt_err_size
  implicit none

  real(dp) :: yields(3), bad_sensitivity(2), returns(3)
  integer :: status
  character(len=160) :: message

  yields = [0.03_dp, 0.031_dp, 0.032_dp]
  bad_sensitivity = 1.0_dp

  returns = total_return(yields, 0.0_dp, status=status, message=message)
  call assert_true(status == tt_err_maturity, 'invalid maturity status')
  call assert_true(all(ieee_is_nan(returns)), 'invalid maturity output')

  returns = total_return(yields, 5.0_dp, scale=0.0_dp, status=status, message=message)
  call assert_true(status == tt_err_scale, 'invalid scale status')

  returns = total_return(yields, 5.0_dp, mdur=bad_sensitivity, &
    status=status, message=message)
  call assert_true(status == tt_err_size, 'invalid sensitivity size status')

  call assert_true(ieee_is_nan(mod_duration(-2.0_dp, 5.0_dp)), &
    'invalid yield duration')
  call assert_true(ieee_is_nan(convexity(0.04_dp, -1.0_dp)), &
    'invalid maturity convexity')

  print '(a)', 'test_validation: PASS'

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      print '(a)', trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_validation
