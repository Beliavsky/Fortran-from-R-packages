! SPDX-License-Identifier: GPL-3.0-only
program test_esreg_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use esback
  implicit none
  real(dp) :: z, loss
  real(dp) :: r(4), q(4), e(4)

  z = -2.0_dp
  call assert_close(g2_curly(z,1), -log(2.0_dp), 1.0e-14_dp, 'G2 curly')
  call assert_close(g2_value(z,1), 0.5_dp, 1.0e-14_dp, 'G2')
  call assert_close(g2_prime(z,1), 0.25_dp, 1.0e-14_dp, 'G2 prime')
  call assert_close(g2_second(z,1), 0.125_dp, 1.0e-14_dp, 'G2 second source convention')

  r = [-2.5_dp, -1.0_dp, -3.0_dp, 0.2_dp]
  q = [-2.0_dp, -1.2_dp, -2.8_dp, -0.5_dp]
  e = [-2.7_dp, -1.8_dp, -3.4_dp, -1.1_dp]
  loss = esr_loss(r, q, e, 0.25_dp)
  call assert_close(loss, 0.6404102953610527_dp, 2.0e-14_dp, 'FZ loss')
  call assert_true(ieee_is_finite(loss), 'finite loss')

  print '(a)', 'test_esreg_core: PASS'
contains
  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(*), intent(in) :: label
    if (abs(actual-expected) > tol) then
      print *, trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not. condition) then
      print *, trim(label)
      error stop 1
    end if
  end subroutine assert_true
end program test_esreg_core
