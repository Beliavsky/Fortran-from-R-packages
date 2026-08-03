! SPDX-License-Identifier: GPL-3.0-only
program esback_demo
  use esback
  implicit none
  integer, parameter :: n = 250
  real(dp) :: r(n), q(n), e(n), s(n), x
  integer :: i
  type(er_backtest_result) :: er
  type(cc_backtest_result) :: cc
  type(esr_backtest_result) :: esr
  type(esreg_options) :: opt

  do i = 1, n
    x = sin(0.071_dp*i) + 0.4_dp*cos(0.19_dp*i)
    s(i) = 0.9_dp + 0.1_dp*sin(0.023_dp*i)
    q(i) = -1.65_dp*s(i) + 0.10_dp*x
    e(i) = -2.05_dp*s(i) + 0.10_dp*x
    r(i) = 0.10_dp*x + s(i)*(0.65_dp*sin(1.71_dp*i) + 0.35_dp*cos(0.37_dp*i))
  end do

  call er_backtest(r, q, e, er, s, 500)
  call cc_backtest(r, q, e, 0.05_dp, cc, s)
  opt%multistarts = 3
  call esr_backtest(r, e, 0.05_dp, 1, esr, b=0, options=opt)

  print '(a,2f12.6)', 'ER p-values: ', er%pvalue_twosided_simple, &
    er%pvalue_twosided_standardized
  print '(a,2f12.6)', 'CC p-values: ', cc%pvalue_twosided_simple, &
    cc%pvalue_twosided_general
  print '(a,f12.6)', 'Strict ESR p-value: ', esr%pvalue_twosided_asymptotic
end program esback_demo
