! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program test_advanced_algorithms
  use kinds_mod, only: dp
  use test_support_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use comoments_mod, only: covariance_matrix, coskewness_matrix, cokurtosis_matrix
  use advanced_moments_mod
  use tail_models_mod
  use analytics_extensions_mod
  implicit none
  integer, parameter :: n=160, p=3
  real(dp) :: r(n,p), f(n,1), rf(n), market(n), asset(n), zinfo(n,1)
  real(dp) :: cov(p,p), target2(p,p), sh2(p,p), sh3(p,p*p), sh4(p,p*p*p)
  real(dp) :: m3(p,p*p), m4(p,p*p*p), ew3(p,p*p), ew4(p,p*p*p)
  real(dp) :: lambda2,lambda3,lambda4,norms,up,down,cap,multi_lambda(2)
  real(dp) :: mt2(p,p),mt3(p,p*p),mt4(p,p*p*p)
  character(len=24) :: target_kinds(2)
  type(mca_result) :: mca3,mca4
  type(nce_result) :: nce,nce_weighted
  type(gpd_fit_result) :: gpd
  type(performance_summary) :: summary
  type(conditional_capm_result) :: ccapm
  real(dp),allocatable :: alpha(:),beta(:),r2(:), losses(:), returns_tail(:)
  real(dp) :: mcvar(p),mces(p),weights(p),vcomp(p),ecomp(p),pvar,pes
  real(dp) :: vse,ese,sse,u,xi,sigma,q,excess,lv,le,kvar,kes,kvc(p),kec(p)
  logical :: ok
  integer :: i,nout,ntail,wins(3),loss_count(3),totals(3),periods(3)
  real(dp) :: outprob(3)

  do i=1,n
    f(i,1)=sin(0.11_dp*real(i,dp))+0.35_dp*cos(0.037_dp*real(i,dp))
    r(i,1)=1.20_dp*f(i,1)+0.16_dp*sin(0.41_dp*real(i,dp))
    r(i,2)=-0.75_dp*f(i,1)+0.20_dp*cos(0.29_dp*real(i,dp))
    r(i,3)=0.45_dp*f(i,1)+0.12_dp*sin(0.23_dp*real(i,dp)+0.4_dp)
  end do
  call covariance_matrix(r,cov,.false.)
  call structured_covariance(r,'indep',target2)
  call assert_close(maxval(abs(target2-diagonal_only(target2))),0.0_dp,1.0e-12_dp,'independent covariance')
  call structured_covariance(r,'observedfactor',target2,f)
  call assert_true(all([(target2(i,i)>0.0_dp,i=1,p)]),'factor covariance positive diagonal')

  call shrink_covariance(r,'constant_correlation',sh2,lambda2)
  call shrink_coskewness(r,'central_symmetric',sh3,lambda3)
  call shrink_cokurtosis(r,'independent',sh4,lambda4)
  call assert_true(lambda2>=0.0_dp .and. lambda2<=1.0_dp,'M2 shrinkage intensity')
  call assert_true(lambda3>=0.0_dp .and. lambda3<=1.0_dp,'M3 shrinkage intensity')
  call assert_true(lambda4>=0.0_dp .and. lambda4<=1.0_dp,'M4 shrinkage intensity')
  call assert_matrix_close(sh2,transpose(sh2),1.0e-10_dp,'shrunk covariance symmetry')
  target_kinds=['independent             ','constant_correlation    ']
  call multi_target_shrink_covariance(r,target_kinds,mt2,multi_lambda)
  call assert_true(all(multi_lambda>=0.0_dp) .and. sum(multi_lambda)<=1.0_dp+1.0e-10_dp,'multi-target M2 weights')
  call multi_target_shrink_coskewness(r,['central_symmetric       ','simaan                  '],mt3,multi_lambda)
  call assert_true(all(multi_lambda>=0.0_dp) .and. sum(multi_lambda)<=1.0_dp+1.0e-10_dp,'multi-target M3 weights')
  call multi_target_shrink_cokurtosis(r,target_kinds,mt4,multi_lambda)
  call assert_true(all(multi_lambda>=0.0_dp) .and. sum(multi_lambda)<=1.0_dp+1.0e-10_dp,'multi-target M4 weights')
  call structured_cokurtosis(r,'constant_correlation',mt4)
  call assert_true(all(ieee_is_finite(mt4)),'constant-correlation M4 target')

  call coskewness_matrix(r,m3);call cokurtosis_matrix(r,m4)
  call ewma_coskewness(r,0.94_dp,ew3);call ewma_cokurtosis(r,0.94_dp,ew4)
  call assert_true(all(ieee_is_finite(ew3)) .and. all(ieee_is_finite(ew4)),'EWMA higher moments finite')

  call m3_mca(r,1,mca3,max_iterations=100,tolerance=1.0e-7_dp)
  call m4_mca(r,1,mca4,max_iterations=100,tolerance=1.0e-7_dp)
  call assert_close(sum(mca3%directions(:,1)**2),1.0_dp,1.0e-7_dp,'M3 MCA direction norm')
  call assert_close(sum(mca4%directions(:,1)**2),1.0_dp,1.0e-7_dp,'M4 MCA direction norm')
  norms=sum(mca3%m3*mca3%m3);call assert_true(norms>0.0_dp,'M3 MCA nonzero reconstruction')
  norms=sum(mca4%m4*mca4%m4);call assert_true(norms>0.0_dp,'M4 MCA nonzero reconstruction')

  call nearest_comoment_estimator(r,1,nce,max_iterations=50,tolerance=1.0e-5_dp)
  call assert_true(nce%objective<huge(1.0_dp)/100.0_dp,'NCE finite objective')
  call assert_true(all(nce%residual_variance>0.0_dp),'NCE positive residual variance')
  call assert_true(all(nce%factor_kurtosis>=nce%factor_skewness**2+1.0_dp),'NCE factor moment inequalities')
  call assert_true(all([(nce%residual_fourth(i)>=nce%residual_variance(i)**2+ &
    nce%residual_third(i)**2/nce%residual_variance(i),i=1,p)]),'NCE residual moment inequalities')
  call assert_matrix_close(nce%covariance,transpose(nce%covariance),1.0e-10_dp,'NCE covariance symmetry')
  call nearest_comoment_estimator(r,1,nce_weighted,max_iterations=20,tolerance=1.0e-4_dp, &
    weight_mode='ridge_diagonal',ridge_alpha=0.2_dp)
  call assert_true(nce_weighted%objective<huge(1.0_dp)/100.0_dp,'weighted NCE finite objective')

  ntail=1000;allocate(losses(ntail),returns_tail(ntail));u=1.0_dp;sigma=0.55_dp;xi=0.20_dp
  do i=1,900
    losses(i)=0.8_dp*real(i-1,dp)/899.0_dp
  end do
  do i=1,100
    q=(real(i,dp)-0.5_dp)/100.0_dp
    excess=sigma/xi*((1.0_dp-q)**(-xi)-1.0_dp)
    losses(900+i)=u+excess
  end do
  returns_tail=-losses/100.0_dp
  call gpd_fit(returns_tail,0.99_dp,gpd,threshold=0.01_dp)
  call assert_true(gpd%converged,'GPD fit convergence')
  call assert_close(gpd%shape,xi,0.18_dp,'GPD shape recovery')
  call assert_close(gpd%scale,sigma/100.0_dp,0.20_dp,'GPD scale recovery')
  call assert_true(gpd%es_value>=gpd%var_value .and. gpd%var_value>gpd%threshold,'GPD risk ordering')
  call assert_true(gpd%var_lower<=gpd%var_value .and. gpd%var_upper>=gpd%var_value,'GPD VaR interval')

  weights=[0.4_dp,0.35_dp,0.25_dp]
  call monte_carlo_asset_risk(r/100.0_dp,0.95_dp,12000,12345_8,mcvar,mces,ok)
  call assert_true(ok .and. all(mcvar>0.0_dp) .and. all(mces>=mcvar),'Monte Carlo asset risk')
  call monte_carlo_portfolio_risk(r/100.0_dp,weights,0.95_dp,16000,54321_8,pvar,pes,vcomp,ecomp,ok)
  call assert_true(ok .and. pvar>0.0_dp .and. pes>=pvar,'Monte Carlo portfolio risk')
  call assert_close(sum(ecomp),pes,0.08_dp,'Monte Carlo ES contributions sum')
  call bootstrap_risk_standard_errors(r(:,1)/100.0_dp,0.95_dp,80,8,9876_8,vse,ese,sse)
  call assert_true(vse>=0.0_dp .and. ese>=0.0_dp .and. sse>=0.0_dp,'bootstrap risk standard errors')
  lv=lognormal_var(r(:,1)/100.0_dp,0.95_dp);le=lognormal_es(r(:,1)/100.0_dp,0.95_dp)
  call assert_true(lv>0.0_dp .and. le>=lv,'lognormal VaR and ES')
  call kernel_portfolio_risk(r/100.0_dp,weights,0.95_dp,kvar,kes,kvc,kec)
  call assert_close(sum(kvc),kvar,1.0e-8_dp,'kernel VaR contributions sum')
  call assert_close(sum(kec),kes,1.0e-8_dp,'kernel ES contributions sum')

  rf=0.0001_dp;market=0.002_dp+0.01_dp*f(:,1);asset=0.0005_dp+1.4_dp*(market-rf)+rf
  call rolling_sfm(asset,market,rf,60,alpha,beta,r2,nout)
  call assert_true(nout==n-59,'rolling CAPM output size')
  call assert_close(beta(nout),1.4_dp,1.0e-7_dp,'rolling CAPM beta')
  call expanding_sfm(asset,market,rf,50,alpha,beta,r2,nout)
  call assert_close(beta(nout),1.4_dp,1.0e-7_dp,'expanding CAPM beta')
  zinfo(:,1)=cos(0.07_dp*[(real(i,dp),i=1,n)])
  asset(1)=rf(1)
  do i=2,n
    asset(i)=rf(i)+0.001_dp+0.4_dp*zinfo(i-1,1)+(0.8_dp+0.3_dp*zinfo(i-1,1))*(market(i)-rf(i))
  end do
  call conditional_capm_fit(asset,market,rf,zinfo,1,ccapm)
  call assert_true(ccapm%success,'conditional CAPM fit')
  call assert_vector_close(ccapm%coefficients,[0.001_dp+0.4_dp*sum(zinfo(:,1))/real(n,dp), &
    0.8_dp+0.3_dp*sum(zinfo(:,1))/real(n,dp),0.4_dp,0.3_dp],1.0e-7_dp,'conditional CAPM coefficients')
  periods=[1,3,12]
  call outperformance_probabilities(r(:,1)/100.0_dp,r(:,2)/100.0_dp,periods,wins,loss_count,totals,outprob)
  call assert_true(all(totals>0) .and. all(outprob>=0.0_dp) .and. all(outprob<=1.0_dp),'outperformance probabilities')
  call capture_ratios(asset,market,up,down,cap)
  call assert_true(ieee_is_finite(up) .and. ieee_is_finite(down) .and. ieee_is_finite(cap),'capture ratios finite')
  call compute_performance_summary(r(:,1)/100.0_dp,12.0_dp,0.95_dp,summary,benchmark=r(:,2)/100.0_dp,rf=0.0_dp,method='historical')
  call assert_true(summary%observations==n .and. summary%var_value>0.0_dp,'performance summary')

  write(*,'(a)')'Advanced moments, NCE/MCA, GPD, Monte Carlo, and dynamic analytics tests passed.'
contains
  pure function diagonal_only(a) result(d)
    real(dp),intent(in)::a(:,:)
    real(dp)::d(size(a,1),size(a,2))
    integer::j
    d=0.0_dp;do j=1,min(size(a,1),size(a,2));d(j,j)=a(j,j);end do
  end function diagonal_only
end program test_advanced_algorithms
