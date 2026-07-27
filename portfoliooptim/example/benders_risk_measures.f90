! SPDX-License-Identifier: GPL-3.0-only
program benders_risk_measures
  use portfoliooptim, only : dp, portfolio_result, bdportfolio_optim, &
    risk_cvar, risk_mad
  implicit none
  real(dp) :: returns(5, 2), probabilities(5), a(2, 2), b(2)
  real(dp) :: lower(2), upper(2)
  type(portfolio_result) :: result
  integer :: code
  character(len=5), parameter :: names(4) = ['CVAR ', 'DCVAR', 'LSAD ', 'MAD  ']

  returns(1, :) = [-0.10_dp, 0.02_dp]
  returns(2, :) = [0.04_dp, -0.03_dp]
  returns(3, :) = [0.08_dp, 0.05_dp]
  returns(4, :) = [0.02_dp, 0.01_dp]
  returns(5, :) = [0.06_dp, 0.03_dp]
  probabilities = 0.20_dp
  a(1, :) = 1.0_dp
  a(2, :) = -1.0_dp
  b = [1.0_dp, -1.0_dp]
  lower = 0.0_dp
  upper = 1.0_dp

  do code = risk_cvar, risk_mad
    result = bdportfolio_optim(returns, probabilities, 0.015_dp, code, &
      0.80_dp, a, b, lower, upper)
    print '(a5,2x,2(f10.6,1x),a,f10.6)', names(code), result%theta, &
      'risk=', result%risk
  end do
end program benders_risk_measures
