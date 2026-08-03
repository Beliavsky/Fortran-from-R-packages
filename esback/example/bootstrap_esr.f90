! SPDX-License-Identifier: GPL-3.0-only
program bootstrap_esr
  use esback
  implicit none
  integer, parameter :: n = 100
  real(dp) :: r(n), e(n)
  integer :: i
  type(esr_backtest_result) :: result
  type(esreg_options) :: opt

  do i = 1, n
    e(i) = -2.0_dp + 0.2_dp*sin(0.08_dp*i)
    r(i) = e(i) + 0.8_dp*sin(0.91_dp*i) + 0.25_dp*cos(1.7_dp*i)
  end do
  opt%multistarts = 1
  opt%sigma_est = sigma_ind
  opt%max_iterations = 1200
  call esr_backtest(r, e, 0.10_dp, 3, result, b=20, options=opt)
  print *, result%pvalue_twosided_asymptotic
  print *, result%pvalue_twosided_bootstrap
  print *, result%pvalue_onesided_bootstrap
end program bootstrap_esr
