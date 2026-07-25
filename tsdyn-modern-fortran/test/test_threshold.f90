! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
program test_threshold
  use tsdyn
  use test_support
  implicit none
  integer,parameter::n=520,k=2
  integer::info,i
  real(dp),allocatable::innov(:,:),y(:,:),y3(:,:),fc(:,:),response(:,:),irf(:,:,:),sresp(:),sirf(:)
  integer,allocatable::reg(:),freg(:)
  real(dp)::stat,thr,spread
  type(tvar_model)::tm,tf,tf3
  type(tvecm_model)::cm,cm3,cmgrid
  type(setar_model)::sm
  real(dp)::beta(2)

  call seed_random(9917)
  allocate(innov(n,k));do i=1,n;innov(i,1)=0.12_dp*random_normal();innov(i,2)=0.10_dp*random_normal();end do

  tm%nvar=2;tm%order=1;tm%nthresh=1;tm%nregime=2;tm%include=include_const;tm%th_delay=1;tm%transition='TAR'
  allocate(tm%thresholds(1),tm%transition_weights(2),tm%coefficients(3,2,2),tm%sigma(2,2))
  tm%thresholds=0.0_dp;tm%transition_weights=[1.0_dp,0.0_dp];tm%coefficients=0.0_dp
  tm%coefficients(1,1,:)=[0.35_dp,0.10_dp];tm%coefficients(1,2,:)=[-0.35_dp,-0.10_dp]
  tm%coefficients(2:3,1,:)=reshape([0.25_dp,0.02_dp,-0.04_dp,0.20_dp],[2,2])
  tm%coefficients(2:3,2,:)=reshape([0.20_dp,-0.03_dp,0.05_dp,0.25_dp],[2,2])
  tm%sigma=0.0_dp;tm%sigma(1,1)=0.0144_dp;tm%sigma(2,2)=0.0100_dp
  call simulate_tvar(tm,n,y,info,innov=innov,start=reshape([-0.2_dp,0.0_dp],[1,2]));call assert_true(info==0,'simulate_tvar')
  call fit_tvar(y,1,include_const,1,tf,info,thresholds=[0.0_dp],weights=[1.0_dp,0.0_dp],trim=0.05_dp)
  call assert_true(info==0,'fit_tvar');call assert_true(all(tf%regime_counts>40),'TVAR regime counts')
  call forecast_tvar(tf,y,5,fc,freg,info);call assert_true(info==0,'forecast_tvar');call assert_true(size(freg)==5,'TVAR forecast regimes')
  call tvar_regimes(tf,y,reg,info);call assert_true(info==0.and.size(reg)==n-1,'tvar_regimes')
  call tvar_lr_statistic(y,1,include_const,stat,thr,info,ngrid=12);call assert_true(info==0,'TVAR LR statistic')
  call regime_irf_tvar(tf,6,1,irf,info,orthogonal=.true.);call assert_true(info==0,'regime TVAR IRF');call assert_all_finite(irf,'TVAR IRF finite')
  call girf_tvar(tf,y(n-2:n,:),[0.2_dp,0.0_dp],5,20,response,info,seed=12);call assert_true(info==0,'TVAR GIRF')
  allocate(y3(n,2));do i=1,n;y3(i,1)=sin(0.05_dp*real(i,dp))+0.05_dp*innov(i,1);y3(i,2)=cos(0.04_dp*real(i,dp))+0.05_dp*innov(i,2);end do
  call fit_tvar(y3,1,include_const,2,tf3,info,thresholds=[-0.30_dp,0.30_dp],weights=[1.0_dp,0.0_dp],trim=0.05_dp)
  call assert_true(info==0,'three-regime TVAR');call assert_true(all(tf3%regime_counts>20),'three-regime TVAR counts')

  ! Cointegrated data with threshold-dependent spread adjustment.
  y=0.0_dp;spread=-0.2_dp
  do i=2,n
    y(i,2)=y(i-1,2)+innov(i,2)
    if(spread<=0.0_dp)then
      spread=0.45_dp*spread+0.08_dp+innov(i,1)
    else
      spread=0.45_dp*spread-0.08_dp+innov(i,1)
    end if
    y(i,1)=y(i,2)+spread
  end do
  beta=[1.0_dp,-1.0_dp]
  call fit_tvecm(y,1,1,include_const,cm,info,beta_fixed=beta,thresholds=[0.0_dp],common='only_ECT',trim_fraction=0.05_dp)
  call assert_true(info==0,'fit_tvecm');call assert_true(all(cm%regime_counts>40),'TVECM regime counts')
  call fit_tvecm(y,1,2,include_const,cm3,info,beta_fixed=beta,thresholds=[-0.08_dp,0.08_dp],common='only_ECT',trim_fraction=0.03_dp)
  call assert_true(info==0,'three-regime TVECM');call assert_true(all(cm3%regime_counts>10),'three-regime TVECM counts')
  call fit_tvecm(y,1,1,include_const,cmgrid,info,thresholds=[0.0_dp],common='only_ECT',trim_fraction=0.05_dp,ngrid_beta=7)
  call assert_true(info==0,'grid-beta TVECM')
  call forecast_tvecm(cm,y,5,fc,freg,info);call assert_true(info==0,'forecast_tvecm')
  call tvecm_regimes(cm,y,reg,info);call assert_true(info==0.and.size(reg)==n-2,'tvecm_regimes')
  call tvecm_lr_statistic(y,1,include_const,beta,stat,thr,info,ngrid=12);call assert_true(info==0,'TVECM LR statistic')
  call girf_tvecm(cm,y(n-3:n,:),[0.2_dp,0.0_dp],5,20,response,info,seed=31);call assert_true(info==0,'TVECM GIRF')

  sm%nregime=2;sm%nthresh=1;sm%pmax=1;sm%include=include_const;sm%th_delay=0;sm%transition='TAR';sm%sigma2=0.04_dp
  allocate(sm%orders(2),sm%thresholds(1),sm%coefficients(2,2));sm%orders=[1,1];sm%thresholds=0.0_dp
  sm%coefficients(:,1)=[0.3_dp,0.25_dp];sm%coefficients(:,2)=[-0.3_dp,0.25_dp]
  call regime_irf_setar(sm,6,1,sirf,info,cumulative=.false.);call assert_true(info==0,'regime SETAR IRF')
  call girf_setar(sm,[-0.2_dp],0.25_dp,5,20,sresp,info,seed=17);call assert_true(info==0,'SETAR GIRF')

  write(*,'(a)') 'Threshold multivariate and nonlinear impulse-response tests passed.'
end program test_threshold
