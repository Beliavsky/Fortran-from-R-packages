! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program example_momentum_backtest
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:),weights(:,:)
  type(backtest_result)::bt
  integer::status
  call generate_sample_prices(156,12,prices,2025_i8)
  call momentum_top_n_strategy(prices,[12.0_dp,5.0_dp,1.0_dp],weights,status)
  if(status/=0)error stop 'strategy failed'
  call run_backtest(prices,weights,100000.0_dp,bt,cost_bps=5.0_dp,frequency=52.0_dp)
  print '(a,f10.4)','total return: ',bt%total_return
  print '(a,f10.4)','annualized Sharpe: ',bt%sharpe
  print '(a,f10.4)','maximum drawdown: ',bt%max_drawdown
end program example_momentum_backtest
