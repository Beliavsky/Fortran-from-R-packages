! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program demo_portfolio_tester
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:),mom(:,:),selection(:,:),weights(:,:),scores(:,:)
  real(dp),allocatable::labels(:,:),features(:,:,:),lag_mom(:,:)
  type(backtest_result)::bt
  type(performance_result)::perf
  integer::status
  call generate_sample_prices(208,15,prices,777_i8)
  call calc_momentum(prices,12,mom)
  call filter_top_n(mom,5,selection)
  call weight_by_rank(selection,mom,weights)
  call run_backtest(prices,weights,100000.0_dp,bt,cost_bps=5.0_dp,frequency=52.0_dp)
  call analyze_performance(bt,52.0_dp,perf)
  call panel_lag(mom,1,lag_mom);call make_labels(prices,4,1,labels)
  allocate(features(208,15,1));features(:,:,1)=lag_mom
  call rolling_fit_predict(features,labels,104,4,4,scores,lambda=0.5_dp)
  status=count(is_finite(scores))
  print '(a)','PortfolioTesteR modern Fortran demonstration'
  print '(a,f10.4)','strategy total return: ',perf%total_return
  print '(a,f10.4)','annualized Sharpe: ',perf%sharpe
  print '(a,f10.4)','maximum drawdown: ',perf%max_drawdown
  print '(a,i0)','rolling ML forecasts: ',status
end program demo_portfolio_tester
