! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program demo_performanceanalytics
  use kinds_mod, only: dp
  use returns_mod, only: portfolio_result, portfolio_returns, annualized_return
  use risk_mod, only: value_at_risk, expected_shortfall
  use drawdown_mod, only: max_drawdown
  use performance_ratios_mod, only: annualized_sharpe_ratio, sortino_ratio
  implicit none
  real(dp)::r(12,3),w(3),rf(12)
  type(portfolio_result)::p
  integer::i
  do i=1,12
    r(i,1)=0.006_dp+0.018_dp*sin(real(i,dp))
    r(i,2)=0.004_dp+0.012_dp*cos(0.7_dp*real(i,dp))
    r(i,3)=0.003_dp+0.008_dp*sin(0.4_dp*real(i,dp))
  end do
  w=[0.5_dp,0.3_dp,0.2_dp];rf=0.0_dp
  call portfolio_returns(r,w,p,rebalance_every=3,transaction_cost=0.001_dp)
  write(*,'(a,es14.6)')'Annualized return: ',annualized_return(p%returns,12.0_dp)
  write(*,'(a,es14.6)')'Annualized Sharpe: ',annualized_sharpe_ratio(p%returns,rf,12.0_dp,.false.)
  write(*,'(a,es14.6)')'Sortino ratio: ',sortino_ratio(p%returns,0.0_dp)
  write(*,'(a,es14.6)')'Modified VaR: ',value_at_risk(p%returns,0.95_dp,'modified')
  write(*,'(a,es14.6)')'Modified ES: ',expected_shortfall(p%returns,0.95_dp,'modified')
  write(*,'(a,es14.6)')'Maximum drawdown: ',max_drawdown(p%returns)
  write(*,'(a,es14.6)')'Ending wealth: ',p%wealth(size(p%wealth))
end program demo_performanceanalytics
