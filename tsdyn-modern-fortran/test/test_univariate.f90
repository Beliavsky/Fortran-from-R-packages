! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
program test_univariate
  use tsdyn
  use test_support
  implicit none
  integer,parameter::n=500
  integer::info,best,i
  integer,allocatable::best_orders(:)
  real(dp),allocatable::x(:),innov(:),fc(:),scores(:),paths(:,:),fitd(:),pvalx(:)
  integer,allocatable::regs(:),nn(:),origins(:)
  real(dp),allocatable::roll(:,:),act(:,:)
  type(ar_model)::ar, ard, ara
  type(setar_model)::sm, sf, sf3, sbest
  type(lstar_model)::lm, lf, lfauto
  type(llar_result)::lr1,lr2
  type(accuracy_metrics)::acc
  type(bbc_test_result)::bbc
  type(kapshin_result)::ks
  real(dp)::pval,obs,stat,thr

  call seed_random(12345)
  allocate(innov(n));do i=1,n;innov(i)=0.25_dp*random_normal();end do
  call simulate_ar([0.65_dp,-0.15_dp],n,x,intercept=0.2_dp,innov=innov,info=info)
  call assert_true(info==0,'simulate_ar')
  call fit_ar(x,2,include_const,'level',ar,info)
  call assert_true(info==0,'fit_ar')
  call assert_close(ar%coefficients(2),0.65_dp,0.12_dp,'AR phi1 recovery')
  call assert_close(ar%coefficients(3),-0.15_dp,0.12_dp,'AR phi2 recovery')
  call forecast_ar(ar,x,5,fc,info);call assert_true(info==0.and.size(fc)==5,'forecast_ar')
  call fit_ar(x,1,include_const,'diff',ard,info);call assert_true(info==0,'AR differences fit')
  call fit_ar(x,1,include_const,'ADF',ara,info);call assert_true(info==0,'AR ADF fit')
  call select_ar_order(x,5,include_const,'BIC',best,scores,info);call assert_true(info==0.and.best>=1.and.best<=3,'select_ar_order')
  acc=compute_accuracy(ar%fitted,x(3:));call assert_finite(acc%rmse,'accuracy')
  call residual_bootstrap_ar(ar,x,4,20,paths,info);call assert_true(info==0.and.all(shape(paths)==[4,20]),'AR bootstrap')

  sm%nregime=2;sm%nthresh=1;sm%pmax=1;sm%include=include_const;sm%th_delay=0;sm%transition='TAR'
  allocate(sm%orders(2),sm%thresholds(1),sm%coefficients(2,2));sm%orders=[1,1];sm%thresholds=0.0_dp
  sm%coefficients(:,1)=[0.5_dp,0.20_dp];sm%coefficients(:,2)=[-0.5_dp,0.20_dp]
  call simulate_setar(sm,n,x,innov=innov,start=[-0.5_dp],info=info);call assert_true(info==0,'simulate_setar')
  call fit_setar(x,[1,1],include_const,1,sf,info,thresholds=[0.0_dp],trim=0.05_dp)
  call assert_true(info==0,'fit_setar fixed threshold')
  call assert_true(all(sf%regime_counts>20),'SETAR regime counts')
  call forecast_setar(sf,x,5,fc,regs,info);call assert_true(info==0,'forecast_setar');call assert_true(size(regs)==5,'forecast SETAR regimes')
  call setar_regimes(sf,x,regs,info);call assert_true(info==0,'setar_regimes');call assert_true(size(regs)==n-1,'setar regime path size')
  call setar_lr_statistic(x,1,include_const,stat,thr,info,ngrid=12);call assert_true(info==0,'setar LR')
  call select_setar_orders(x,2,include_const,'BIC',best_orders,sbest,info,ngrid=10);call assert_true(info==0.and.size(best_orders)==2,'SETAR order selection')

  lm%p_low=1;lm%p_high=1;lm%pmax=1;lm%include=include_const;lm%th_delay=0;lm%gamma=5.0_dp;lm%threshold=0.0_dp
  allocate(lm%coef_low(2),lm%coef_high(2));lm%coef_low=[-0.15_dp,0.6_dp];lm%coef_high=[0.30_dp,-0.25_dp]
  call simulate_lstar(lm,n,x,innov=innov,start=[0.0_dp],info=info);call assert_true(info==0,'simulate_lstar')
  call fit_lstar(x,1,1,include_const,lf,info,th_delay=0,gamma=5.0_dp,threshold=0.0_dp)
  call assert_true(info==0,'fit_lstar');call assert_true(lf%gamma>0.0_dp,'LSTAR gamma')
  call fit_lstar(x,1,1,include_const,lfauto,info,th_delay=0,ngamma=5,nthreshold=10)
  call assert_true(info==0.and.lfauto%gamma>0.0_dp,'automatic LSTAR search')
  call forecast_lstar(lf,x,5,fc,info);call assert_true(info==0,'forecast_lstar')
  call assert_true(lstar_transition(-1.0_dp,5.0_dp,0.0_dp)<lstar_transition(1.0_dp,5.0_dp,0.0_dp),'LSTAR transition')

  do i=1,n;x(i)=sin(0.08_dp*real(i,dp))+0.03_dp*innov(i);end do
  call fit_setar(x,[1,1,1],include_const,2,sf3,info,thresholds=[-0.30_dp,0.30_dp],trim=0.05_dp)
  call assert_true(info==0,'three-regime SETAR');call assert_true(all(sf3%regime_counts>20),'three-regime SETAR counts')
  call llar_fit_curve(x,2,1,2,[0.12_dp,0.20_dp],lr1,info,'direct');call assert_true(info==0,'LLAR direct')
  call llar_fit_curve(x,2,1,2,[0.12_dp,0.20_dp],lr2,info,'box');call assert_true(info==0,'LLAR box')
  call assert_true(all(lr1%usable_points==lr2%usable_points),'LLAR direct/box usable points')
  call assert_true(maxval(abs(lr1%normalized_error-lr2%normalized_error))<1.0e-10_dp,'LLAR direct/box equivalence')
  call llar_fitted(x,2,1,2,0.20_dp,fitd,nn,info,'auto');call assert_true(info==0,'LLAR fitted');call assert_true(count(nn>0)>100,'LLAR neighbors')
  call llar_predict(x,2,1,2,0.20_dp,3,fc,info,'auto');call assert_true(info==0.and.size(fc)==3,'LLAR predict')

  call assert_true(correlation_integral(x,2,1,0.3_dp,2)>=0.0_dp,'correlation integral')
  call assert_finite(delta_statistic(x,2,1,0.3_dp),'delta statistic')
  call assert_finite(delta_linear_statistic(x,2,1),'delta linear statistic')
  call delta_shuffle_test(x,2,1,0.3_dp,9,pval,obs,info,seed=44);call assert_true(info==0,'delta shuffle');call assert_true(pval>0.0_dp.and.pval<=1.0_dp,'delta shuffle p-value')
  call delta_linearity_test(x,2,1,0.3_dp,5,pval,obs,info,seed=44,pmax=3);call assert_true(info==0,'delta linearity')
  call bbc_unit_root_test(x,1,0.15_dp,bbc,info,ngrid=8);call assert_true(info==0,'BBC test')
  call kapshin_test(x,1,include_const,1.0_dp,0.5_dp,10,ks,info,npoints=7);call assert_true(info==0,'KapShin test');call assert_true(ks%valid_pairs>0,'KapShin pairs')

  call rolling_forecast_ar(x,350,2,2,include_const,roll,act,origins,info,window=250)
  call assert_true(info==0,'rolling AR');call assert_true(size(roll,1)==2.and.size(roll,2)==149,'rolling AR shape')
  call rolling_forecast_setar(x,400,1,[1,1],include_const,1,roll,act,origins,info,thresholds=[0.0_dp],window=300)
  call assert_true(info==0,'rolling SETAR');call assert_true(size(roll,2)==100,'rolling SETAR shape')
  call resample_vector(x,'block',pvalx,info,block_length=12);call assert_true(info==0,'block resample');call assert_true(size(pvalx)==n,'block resample size')
  call resample_vector(x,'wild2',pvalx,info);call assert_true(info==0,'wild resample')

  write(*,'(a)') 'Univariate and nonlinear tests passed.'
end program test_univariate
