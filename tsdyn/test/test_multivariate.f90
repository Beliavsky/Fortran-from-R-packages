! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
program test_multivariate
  use tsdyn
  use test_support
  implicit none
  integer,parameter::n=520,k=2
  integer::info,best,i
  real(dp),allocatable::innov(:,:),y(:,:),fc(:,:),irf(:,:,:),fevd(:,:,:),scores(:),paths(:,:,:)
  real(dp),allocatable::a(:,:,:),trace(:),maxeig(:),rolling(:,:,:),actual(:,:,:),bres(:,:)
  integer,allocatable::origins(:)
  real(dp)::spread
  type(var_model)::vm,vmd,vma
  type(vecm_model)::cm,cmml,cm2ols
  real(dp)::coef(3,2),beta(2,1)

  call seed_random(883)
  allocate(innov(n,k));do i=1,n;innov(i,1)=0.16_dp*random_normal();innov(i,2)=0.13_dp*random_normal();end do
  coef(1,:)=[0.10_dp,-0.05_dp]
  coef(2,:)=[0.55_dp,0.10_dp]
  coef(3,:)=[-0.08_dp,0.42_dp]
  call simulate_var(coef,1,include_const,n,y,info,innov=innov,start=reshape([0.0_dp,0.0_dp],[1,2]))
  call assert_true(info==0,'simulate_var')
  call fit_var(y,1,include_const,vm,info);call assert_true(info==0,'fit_var')
  call fit_var(y,1,include_const,vmd,info,specification='diff');call assert_true(info==0,'VAR differences fit')
  call fit_var(y,1,include_const,vma,info,specification='ADF');call assert_true(info==0,'VAR ADF fit')
  call assert_close(vm%coefficients(2,1),coef(2,1),0.08_dp,'VAR coefficient recovery')
  call forecast_var(vm,y,6,fc,info);call assert_true(info==0,'forecast_var');call assert_true(all(shape(fc)==[6,2]),'VAR forecast shape')
  call impulse_response_var(vm,8,irf,info,orthogonal=.true.);call assert_true(info==0,'VAR IRF');call assert_all_finite(irf,'VAR IRF finite')
  call fevd_var(vm,8,fevd,info);call assert_true(info==0,'VAR FEVD');call assert_true(maxval(abs(sum(fevd,dim=3)-1.0_dp))<1.0e-10_dp,'FEVD rows sum to one')
  call select_var_lag(y,4,include_const,'BIC',best,scores,info);call assert_true(info==0.and.best>=1.and.best<=2,'VAR lag selection')
  call bootstrap_var(vm,y,5,20,paths,info,wild=.false.);call assert_true(info==0,'VAR residual bootstrap');call assert_true(all(shape(paths)==[5,2,20]),'VAR bootstrap shape')
  call bootstrap_var(vm,y,3,10,paths,info,wild=.true.);call assert_true(info==0,'VAR wild bootstrap')
  call rolling_forecast_var(y,400,2,1,include_const,rolling,actual,origins,info,window=300)
  call assert_true(info==0,'rolling VAR');call assert_true(size(rolling,3)==119,'rolling VAR origins')
  call block_resample_matrix(y,15,bres,info);call assert_true(info==0.and.all(shape(bres)==[n,k]),'matrix block resample')

  ! Construct a cointegrated pair: y1-y2 is stationary while y2 has a common stochastic trend.
  y=0.0_dp;spread=0.0_dp
  do i=2,n
    y(i,2)=y(i-1,2)+innov(i,2)
    spread=0.55_dp*spread+innov(i,1)
    y(i,1)=y(i,2)+spread
  end do
  beta(:,1)=[1.0_dp,-1.0_dp]
  call fit_vecm(y,1,1,include_const,'fixed',cm,info,beta_fixed=beta);call assert_true(info==0,'fixed-beta VECM')
  call fit_vecm(y,1,1,include_const,'2OLS',cm2ols,info);call assert_true(info==0,'two-step VECM')
  call assert_true(cm%alpha(1,1)<0.0_dp,'VECM error correction sign')
  call forecast_vecm(cm,y,5,fc,info);call assert_true(info==0,'forecast_vecm')
  call vecm_var_coefficients(cm,a,info);call assert_true(info==0.and.size(a,3)==2,'VECM to VAR')
  call impulse_response_vecm(cm,6,irf,info,orthogonal=.true.);call assert_true(info==0,'VECM IRF')
  call simulate_vecm(cm,40,bres,info,innov=innov(1:40,:),start=y(1:2,:));call assert_true(info==0,'simulate_vecm')
  call select_vecm_rank(y,1,include_const,1,'BIC',best,scores,info);call assert_true(info==0.and.best==1,'VECM rank selection')

  call fit_vecm(y,1,1,include_const,'ML',cmml,info);call assert_true(info==0,'Johansen ML VECM')
  call johansen_statistics(cmml,trace,maxeig,info);call assert_true(info==0,'Johansen statistics')
  call assert_all_finite(trace,'Johansen trace finite');call assert_all_finite(maxeig,'Johansen max eigen finite')

  write(*,'(a)') 'Multivariate linear and cointegration tests passed.'
end program test_multivariate
