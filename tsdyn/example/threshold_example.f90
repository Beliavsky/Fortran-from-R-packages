! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
program threshold_example
  use tsdyn
  implicit none
  integer,parameter::n=500
  integer::i,info
  real(dp),allocatable::x(:),e(:),fc(:)
  integer,allocatable::reg(:)
  type(setar_model)::truth,fit
  call seed_random(77)
  truth%nregime=2;truth%nthresh=1;truth%pmax=1;truth%include=include_const;truth%th_delay=0;truth%transition='TAR'
  allocate(truth%orders(2),truth%thresholds(1),truth%coefficients(2,2));truth%orders=[1,1];truth%thresholds=0.0_dp
  truth%coefficients(:,1)=[0.45_dp,0.2_dp];truth%coefficients(:,2)=[-0.45_dp,0.2_dp]
  allocate(e(n));do i=1,n;e(i)=0.20_dp*random_normal();end do
  call simulate_setar(truth,n,x,innov=e,start=[-0.2_dp],info=info);if(info/=0)error stop 'simulation failed'
  call fit_setar(x,[1,1],include_const,1,fit,info,thresholds=[0.0_dp],trim=0.05_dp);if(info/=0)error stop 'fit failed'
  call forecast_setar(fit,x,8,fc,reg,info);if(info/=0)error stop 'forecast failed'
  write(*,'(a,2i8)') 'Regime observations:',fit%regime_counts
  write(*,'(a,8f9.4)') 'Forecast:',fc
  write(*,'(a,8i5)') 'Regimes: ',reg
end program threshold_example
