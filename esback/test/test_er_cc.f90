! SPDX-License-Identifier: GPL-3.0-only
program test_er_cc
  use esback
  implicit none
  integer, parameter :: n = 40
  real(dp) :: r(n), q(n), e(n), s(n)
  integer :: i
  type(er_backtest_result) :: er
  type(cc_backtest_result) :: cc

  do i = 1, n
    r(i) = -0.02_dp + 0.012_dp*sin(0.73_dp*i) + 0.006_dp*cos(0.19_dp*i)
    q(i) = -0.021_dp + 0.002_dp*sin(0.17_dp*i)
    e(i) = q(i) - 0.008_dp - 0.001_dp*cos(0.31_dp*i)
    s(i) = 0.012_dp + 0.001_dp*sin(0.11_dp*i)
  end do

  call er_backtest(r, q, e, er, s, 250)
  call assert_true(er%status == esback_ok, 'ER status')
  call assert_true(er%n_exceedances == 16, 'ER exceedance count')
  call assert_close(er%pvalue_twosided_simple, 0.836_dp, 1.0e-14_dp, 'ER simple two-sided')
  call assert_close(er%pvalue_onesided_simple, 0.448_dp, 1.0e-14_dp, 'ER simple one-sided')
  call assert_close(er%pvalue_twosided_standardized, 0.892_dp, 1.0e-14_dp, 'ER standardized two-sided')
  call assert_close(er%pvalue_onesided_standardized, 0.460_dp, 1.0e-14_dp, 'ER standardized one-sided')

  call cc_backtest(r, q, e, 0.1_dp, cc, s, .true.)
  call assert_true(cc%status == esback_ok, 'CC status')
  call assert_close(cc%pvalue_twosided_simple, 0.0042733555437817665_dp, 5.0e-13_dp, 'CC simple two-sided')
  call assert_close(cc%pvalue_onesided_simple, 0.006380448290648433_dp, 5.0e-13_dp, 'CC simple one-sided')
  call assert_close(cc%pvalue_twosided_general, 0.8787136534102065_dp, 5.0e-13_dp, 'CC general two-sided')
  call assert_close(cc%pvalue_onesided_general, 0.008897508299635554_dp, 5.0e-13_dp, 'CC general one-sided')

  print '(a)', 'test_er_cc: PASS'
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
end program test_er_cc
