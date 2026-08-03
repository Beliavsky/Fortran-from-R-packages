! SPDX-License-Identifier: GPL-2.0-or-later
program test_gls
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_MAX_ITER
  use nlme_types
  use nlme_gls, only : fit_gls
  use test_support
  implicit none
  integer,parameter::n=80
  real(dp)::y(n),x(n,2),time(n),e(n),innov
  integer::group(n),i
  type(correlation_spec)::corr
  type(gls_result)::fit
  type(nlme_control)::ctl
  x(:,1)=1.0_dp
  time(1)=1.0_dp;x(1,2)=(1.0_dp-40.5_dp)/40.0_dp
  innov=0.12_dp*sin(1.71_dp)+0.08_dp*cos(0.43_dp);e(1)=innov
  do i=2,n
    time(i)=real(i,dp);x(i,2)=(real(i,dp)-40.5_dp)/40.0_dp
    innov=0.12_dp*sin(1.71_dp*real(i,dp))+0.08_dp*cos(0.43_dp*real(i,dp))
    e(i)=0.65_dp*e(i-1)+innov
  end do
  y=1.5_dp-0.8_dp*x(:,2)+e;group=1
  corr%kind=COR_AR1;allocate(corr%par(1));corr%par=0.2_dp
  ctl%reml=.true.;ctl%max_iter=500;ctl%tolerance=1.0e-8_dp;ctl%step=0.2_dp
  call fit_gls(y,x,fit,corr,time=time,group=group,control=ctl)
  call check(fit%status==NLME_SUCCESS .or. fit%status==NLME_MAX_ITER,'GLS status')
  call check_close(fit%beta(1),1.5_dp,0.08_dp,'GLS intercept')
  call check_close(fit%beta(2),-0.8_dp,0.10_dp,'GLS slope')
  call check(fit%correlation_parameters(1)>0.2_dp .and. fit%correlation_parameters(1)<0.95_dp,'GLS rho')
  call check(fit%sigma>0.0_dp,'GLS sigma')
  call check(size(fit%residuals)==n,'GLS residual size')
  write(*,'(a)')'test_gls: PASS'
end program test_gls
