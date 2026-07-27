! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
program test_metrics
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pbo, only : dp, column_mean, column_sum, sharpe_ratio, omega_ratio
  implicit none
  real(dp) :: x(4,2), values(2)

  x(:,1) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
  x(:,2) = [2.0_dp, 2.0_dp, -1.0_dp, -1.0_dp]
  call column_mean(x,values)
  call assert_close(values(1),0.5_dp,1.0e-14_dp,'column mean 1')
  call assert_close(values(2),0.5_dp,1.0e-14_dp,'column mean 2')
  call column_sum(x,values)
  call assert_close(values(1),2.0_dp,1.0e-14_dp,'column sum 1')
  call assert_close(values(2),2.0_dp,1.0e-14_dp,'column sum 2')
  call sharpe_ratio(x,values)
  call assert_close(values(1),0.3872983346207417_dp,1.0e-14_dp,'sharpe 1')
  call assert_close(values(2),0.2886751345948129_dp,1.0e-14_dp,'sharpe 2')
  call omega_ratio(x,values)
  call assert_close(values(1),3.0_dp,1.0e-14_dp,'omega 1')
  call assert_close(values(2),2.0_dp,1.0e-14_dp,'omega 2')
  call assert_true(all(ieee_is_finite(values)), 'finite metrics')
  print '(a)', 'test_metrics: PASS'
contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual-expected) > tolerance) then
      print '(a,2es24.15)', 'FAILED: '//trim(label)//' actual/expected ',actual,expected
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAILED: '//trim(label)
      error stop 1
    end if
  end subroutine assert_true
end program test_metrics
