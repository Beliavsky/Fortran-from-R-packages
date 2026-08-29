! SPDX-License-Identifier: GPL-2.0-or-later
program test_lme
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_MAX_ITER
  use nlme_types
  use nlme_lme, only : fit_lme
  use test_support
  implicit none
  integer,parameter::ng=12,m=7,n=ng*m
  real(dp)::y(n),x(n,2),z(n,1),time(n),btrue(ng),tt,eps
  integer::group(n),i,g,k
  type(pd_spec)::random
  type(lme_result)::fit
  type(nlme_control)::ctl
  k=0
  do g=1,ng
    btrue(g)=0.45_dp*sin(1.3_dp*real(g,dp))
    do i=1,m
      k=k+1
      tt=real(i-1,dp)/real(m-1,dp)
      x(k,:)=[1.0_dp,tt]
      z(k,1)=1.0_dp
      group(k)=g
      time(k)=real(i,dp)
      eps=0.08_dp*sin(2.1_dp*real(k,dp))+0.03_dp*cos(0.7_dp*real(k,dp))
      y(k)=2.0_dp+1.25_dp*tt+btrue(g)+eps
    end do
  end do
  random%kind=PD_DIAG
  random%dim=1
  allocate(random%par(1))
  random%par=log(0.25_dp)
  ctl%reml=.true.
  ctl%max_iter=700
  ctl%tolerance=2.0e-7_dp
  ctl%step=0.15_dp
  call fit_lme(y,x,z,group,fit,random=random,time=time,control=ctl)
  call check(fit%status==NLME_SUCCESS .or. fit%status==NLME_MAX_ITER,'LME status')
  call check_close(fit%beta(1),2.0_dp,0.10_dp,'LME intercept')
  call check_close(fit%beta(2),1.25_dp,0.08_dp,'LME slope')
  call check(fit%random_covariance(1,1)>0.02_dp,'LME random variance')
  call check(fit%sigma>0.0_dp .and. fit%sigma<0.3_dp,'LME residual sigma')
  call check(sum(fit%random_effects(:,1)*btrue)>0.0_dp,'LME BLUP direction')
  call check(sum(fit%residuals_conditional**2)<sum(fit%residuals_marginal**2),'conditional residual improvement')
  write(*,'(a)')'test_lme: PASS'
end program test_lme
