! SPDX-License-Identifier: GPL-2.0-or-later
program performance_indices
   use jfe
   implicit none

   real(dp), parameter :: returns(8) = [0.012_dp, -0.008_dp, 0.015_dp, 0.004_dp, &
      -0.011_dp, 0.009_dp, 0.006_dp, -0.003_dp]

   print '(a,f12.6)', 'Annualized return: ', return_annualized(returns, 12.0_dp)
   print '(a,f12.6)', 'Annualized Sharpe: ', sharpe_ratio_annualized(returns, scale=12.0_dp)
   print '(a,f12.6)', 'Sortino ratio:     ', sortino_ratio(returns)
   print '(a,f12.6)', 'Maximum drawdown:  ', max_drawdown(returns)
   print '(a,f12.6)', 'Calmar ratio:      ', calmar_ratio(returns, 12.0_dp)
end program performance_indices
