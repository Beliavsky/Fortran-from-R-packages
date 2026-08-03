! SPDX-License-Identifier: GPL-3.0-only
program test_esr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use esback
  implicit none
  integer, parameter :: n = 120
  real(dp) :: r(n), q(n), e(n), x
  integer :: i, version
  type(esr_backtest_result) :: result
  type(esreg_options) :: opt

  do i = 1, n
    x = sin(0.37_dp*i) + 0.5_dp*cos(0.11_dp*i)
    r(i) = 0.7_dp*x - 1.8_dp + 0.25_dp*sin(1.7_dp*i)
    q(i) = 0.7_dp*x - 1.95_dp
    e(i) = 0.7_dp*x - 2.25_dp
  end do

  opt%multistarts = 2
  opt%max_iterations = 1500
  opt%sigma_est = sigma_ind

  call esr_backtest(r, e, 0.1_dp, 1, result, b=0, options=opt)
  call assert_true(result%status == esback_ok, 'strict ESR status')
  call assert_close(result%fit%coefficients_q(1), 0.2279928_dp, 2.0e-3_dp, 'strict q intercept')
  call assert_close(result%fit%coefficients_q(2), 1.0073121_dp, 2.0e-3_dp, 'strict q slope')
  call assert_close(result%fit%coefficients_e(1), 0.2291886_dp, 3.0e-3_dp, 'strict ES intercept')
  call assert_close(result%fit%coefficients_e(2), 1.0108174_dp, 3.0e-3_dp, 'strict ES slope')
  call assert_true(result%fit%covariance_available, 'strict covariance')
  call assert_true(all(ieee_is_finite(result%fit%covariance)), 'finite strict covariance')

  do version = 2, 3
    if (version == 2) then
      call esr_backtest(r, e, 0.1_dp, version, result, q=q, b=0, options=opt)
    else
      call esr_backtest(r, e, 0.1_dp, version, result, b=0, options=opt)
    end if
    call assert_true(result%status == esback_ok, 'ESR version status')
    call assert_true(ieee_is_finite(result%statistic), 'finite ESR statistic')
    call assert_true(result%pvalue_twosided_asymptotic >= 0.0_dp .and. &
      result%pvalue_twosided_asymptotic <= 1.0_dp, 'valid ESR p-value')
  end do

  print '(a)', 'test_esr: PASS'
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
end program test_esr
