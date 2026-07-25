! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
program demo_tsdyn
  use tsdyn
  implicit none
  integer,parameter::n=400
  integer::i,info
  real(dp),allocatable::e(:),x(:),fc(:),y(:,:),mvinnov(:,:),mvfc(:,:),irf(:,:,:)
  type(ar_model)::ar
  type(setar_model)::tar
  type(var_model)::vm
  real(dp)::coef(3,2)
  call seed_random(20260724)
  allocate(e(n));do i=1,n;e(i)=0.25_dp*random_normal();end do
  call simulate_ar([0.6_dp,-0.15_dp],n,x,intercept=0.1_dp,innov=e,info=info)
  if(info/=0)error stop 'AR simulation failed'
  call fit_ar(x,2,include_const,'level',ar,info);if(info/=0)error stop 'AR fit failed'
  call forecast_ar(ar,x,5,fc,info);if(info/=0)error stop 'AR forecast failed'
  write(*,'(a,3f10.4)') 'AR coefficients: ',ar%coefficients
  write(*,'(a,5f10.4)') 'AR forecast:     ',fc

  call fit_setar(x,[1,1],include_const,1,tar,info,ngrid=15,trim=0.10_dp)
  if(info==0)write(*,'(a,f10.4,a,2i7)') 'SETAR threshold: ',tar%thresholds(1),' counts:',tar%regime_counts

  allocate(mvinnov(n,2));do i=1,n;mvinnov(i,1)=0.15_dp*random_normal();mvinnov(i,2)=0.12_dp*random_normal();end do
  coef(1,:)=[0.05_dp,-0.02_dp];coef(2,:)=[0.55_dp,0.08_dp];coef(3,:)=[-0.06_dp,0.40_dp]
  call simulate_var(coef,1,include_const,n,y,info,innov=mvinnov,start=reshape([0.0_dp,0.0_dp],[1,2]));if(info/=0)error stop 'VAR simulation failed'
  call fit_var(y,1,include_const,vm,info);if(info/=0)error stop 'VAR fit failed'
  call forecast_var(vm,y,3,mvfc,info);if(info/=0)error stop 'VAR forecast failed'
  call impulse_response_var(vm,5,irf,info,orthogonal=.true.);if(info/=0)error stop 'VAR IRF failed'
  write(*,'(a,2f10.4)') 'VAR first forecast: ',mvfc(1,:)
  write(*,'(a,2f10.4)') 'VAR impact response:',irf(0,1,:)
end program demo_tsdyn
