! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program portfolio_contributions
  use kinds_mod, only: dp
  use returns_mod, only: portfolio_result, portfolio_returns
  implicit none
  real(dp)::r(4,2),w(2)
  type(portfolio_result)::result
  integer::i
  r=reshape([0.02_dp,-0.01_dp,0.03_dp,0.01_dp,0.01_dp,0.02_dp,-0.02_dp,0.04_dp],[4,2])
  w=[0.6_dp,0.4_dp]
  call portfolio_returns(r,w,result,rebalance_every=2)
  write(*,'(a)')'period return contribution_1 contribution_2 turnover'
  do i=1,4
    write(*,'(i3,4(1x,f12.7))')i,result%returns(i),result%contributions(i,1), &
      result%contributions(i,2),result%turnover(i)
  end do
end program portfolio_contributions
