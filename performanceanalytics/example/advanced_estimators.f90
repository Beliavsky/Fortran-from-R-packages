! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program advanced_estimators
  use kinds_mod, only: dp
  use advanced_moments_mod
  use tail_models_mod
  implicit none
  integer,parameter::n=240,p=3
  real(dp)::r(n,p),weights(p),mcvar(p),mces(p),vc(p),ec(p),pv,pe
  type(nce_result)::nce
  type(mca_result)::mca
  type(gpd_fit_result)::gpd
  logical::ok
  integer::i
  do i=1,n
    r(i,1)=0.007_dp*sin(0.09_dp*real(i,dp))+0.003_dp*cos(0.31_dp*real(i,dp))
    r(i,2)=-0.004_dp*sin(0.09_dp*real(i,dp))+0.005_dp*cos(0.17_dp*real(i,dp))
    r(i,3)=0.003_dp*sin(0.09_dp*real(i,dp))+0.004_dp*sin(0.23_dp*real(i,dp))
  end do
  call m3_mca(r,1,mca)
  call nearest_comoment_estimator(r,1,nce,max_iterations=40)
  call gpd_fit(r(:,1),0.95_dp,gpd,p_threshold=0.90_dp)
  weights=[0.4_dp,0.35_dp,0.25_dp]
  call monte_carlo_asset_risk(r,0.95_dp,10000,20260725_8,mcvar,mces,ok)
  call monte_carlo_portfolio_risk(r,weights,0.95_dp,12000,20260726_8,pv,pe,vc,ec,ok)
  write(*,'(a,l1,a,i0)')'M3 MCA converged: ',mca%converged,', iterations: ',mca%iterations
  write(*,'(a,es12.4)')'NCE objective: ',nce%objective
  write(*,'(a,2(1x,es12.4))')'GPD VaR/ES:',gpd%var_value,gpd%es_value
  write(*,'(a,2(1x,es12.4))')'Monte Carlo portfolio VaR/ES:',pv,pe
end program advanced_estimators
