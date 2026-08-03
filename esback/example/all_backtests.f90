! SPDX-License-Identifier: GPL-3.0-only
program all_backtests
  use esback
  implicit none
  integer, parameter :: n = 160
  real(dp) :: r(n), q(n), e(n), s(n)
  integer :: i, version
  type(er_backtest_result) :: er
  type(cc_backtest_result) :: cc
  type(esr_backtest_result) :: esr
  type(esreg_options) :: opt

  do i = 1, n
    s(i) = 0.8_dp + 0.1_dp*cos(0.04_dp*i)
    q(i) = -1.5_dp*s(i)
    e(i) = -1.9_dp*s(i)
    r(i) = s(i)*(sin(0.83_dp*i) + 0.25_dp*cos(2.1_dp*i))
  end do
  opt%multistarts = 2
  opt%sigma_est = sigma_ind

  call er_backtest(r, q, e, er, s, 250)
  call cc_backtest(r, q, e, 0.10_dp, cc, s)
  print *, 'ER:', er%pvalue_twosided_simple, er%pvalue_twosided_standardized
  print *, 'CC:', cc%pvalue_twosided_simple, cc%pvalue_twosided_general

  do version = 1, 3
    if (version == 2) then
      call esr_backtest(r, e, 0.10_dp, version, esr, q=q, options=opt)
    else
      call esr_backtest(r, e, 0.10_dp, version, esr, options=opt)
    end if
    print *, 'ESR version', version, esr%pvalue_twosided_asymptotic
  end do
end program all_backtests
