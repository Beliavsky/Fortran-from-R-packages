! SPDX-License-Identifier: GPL-3.0-only
program test_benders
  use portfoliooptim, only : dp, portfolio_result, bdportfolio_optim, &
    risk_cvar, risk_dcvar, risk_lsad, risk_mad
  implicit none
  real(dp) :: returns(5, 2), probabilities(5), aconstr(2, 2), bconstr(2)
  real(dp) :: lower(2), upper(2)
  type(portfolio_result) :: result

  returns(1, :) = [-0.10_dp, 0.02_dp]
  returns(2, :) = [0.04_dp, -0.03_dp]
  returns(3, :) = [0.08_dp, 0.05_dp]
  returns(4, :) = [0.02_dp, 0.01_dp]
  returns(5, :) = [0.06_dp, 0.03_dp]
  probabilities = 0.20_dp
  aconstr(1, :) = 1.0_dp
  aconstr(2, :) = -1.0_dp
  bconstr = [1.0_dp, -1.0_dp]
  lower = 0.0_dp
  upper = 1.0_dp

  result = bdportfolio_optim(returns, probabilities, 0.015_dp, risk_cvar, &
    0.80_dp, aconstr, bconstr, lower, upper, 200, 1.0e-10_dp)
  call assert_true(result%converged)
  call assert_close(result%theta(1), 0.2631578947368421_dp, 1.0e-10_dp)
  call assert_close(result%risk, 0.0115789473684211_dp, 1.0e-10_dp)

  result = bdportfolio_optim(returns, probabilities, 0.015_dp, risk_dcvar, &
    0.80_dp, aconstr, bconstr, lower, upper, 200, 1.0e-10_dp)
  call assert_true(result%converged)
  call assert_close(result%risk, 0.0286315789473684_dp, 1.0e-10_dp)

  result = bdportfolio_optim(returns, probabilities, 0.015_dp, risk_lsad, &
    0.80_dp, aconstr, bconstr, lower, upper, 200, 1.0e-10_dp)
  call assert_true(result%converged)
  call assert_close(result%theta(1), 0.0322580645161290_dp, 1.0e-10_dp)
  call assert_close(result%risk, 0.00993548387096774_dp, 1.0e-10_dp)

  result = bdportfolio_optim(returns, probabilities, 0.015_dp, risk_mad, &
    0.80_dp, aconstr, bconstr, lower, upper, 200, 1.0e-10_dp)
  call assert_true(result%converged)
  call assert_close(result%risk, 0.0198709677419355_dp, 1.0e-10_dp)
  call assert_close(sum(result%theta), 1.0_dp, 1.0e-10_dp)
  print '(a)', 'test_benders: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_benders
