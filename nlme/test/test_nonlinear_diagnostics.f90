! SPDX-License-Identifier: GPL-2.0-or-later
program test_nonlinear_diagnostics
  use nlme
  use test_support
  implicit none
  integer, parameter :: n=60, ng=6, m=8, nn=ng*m
  real(dp) :: xd(n,1),y(n),theta0(3)
  real(dp) :: xnm(nn,1),ynm(nn),true_b(ng),tt
  real(dp), allocatable :: acf(:),ysim(:),means(:),sds(:),coefs(:,:),sigmas(:),nlspars(:,:),nlssigma(:)
  integer, allocatable :: pairs(:),levels(:),counts(:)
  integer :: i,g,k,status,groups(nn),ri(1)
  type(nonlinear_result) :: gfit,nfit
  type(pd_spec) :: random
  type(nlme_control) :: ctl
  type(variogram_result) :: vg
  real(dp) :: xs(nn,2),zs(nn,1),gmat(1,1)

  do i=1,n
    xd(i,1)=real(i-1,dp)/20.0_dp
  end do
  y=2.2_dp*exp(-0.7_dp*xd(:,1))+0.4_dp+0.02_dp*sin([(real(i,dp),i=1,n)])
  theta0=[1.5_dp,0.4_dp,0.2_dp]
  ctl%max_iter=200;ctl%tolerance=1.0e-8_dp
  call fit_gnls(exponential_decay_model,y,xd,theta0,gfit,control=ctl)
  call check(gfit%status==NLME_SUCCESS,'GNLS status')
  call check_close(gfit%parameters(1),2.2_dp,0.04_dp,'GNLS amplitude')
  call check_close(gfit%parameters(2),0.7_dp,0.05_dp,'GNLS rate')
  call check_close(gfit%parameters(3),0.4_dp,0.06_dp,'GNLS offset')

  k=0
  do g=1,ng
    true_b(g)=0.25_dp*sin(real(g,dp))
    do i=1,m
      k=k+1;tt=real(i-1,dp)/4.0_dp
      xnm(k,1)=tt;groups(k)=g
      ynm(k)=(1.8_dp+true_b(g))*exp(-0.55_dp*tt)+0.25_dp+0.015_dp*cos(1.7_dp*real(k,dp))
    end do
  end do
  ri=[1]
  random%kind=PD_DIAG;random%dim=1;random%fixed=.true.;allocate(random%par(1));random%par=log(0.22_dp)
  theta0=[1.5_dp,0.4_dp,0.1_dp]
  ctl%max_outer=12;ctl%max_iter=120;ctl%tolerance=2.0e-5_dp
  call fit_nlme(exponential_decay_model,ynm,xnm,groups,theta0,ri,nfit,random=random,control=ctl,fixed_sigma=0.03_dp)
  call check(nfit%status==NLME_SUCCESS .or. nfit%status==NLME_MAX_ITER,'NLME status')
  call check_close(nfit%parameters(1),1.8_dp,0.15_dp,'NLME amplitude')
  call check_close(nfit%parameters(2),0.55_dp,0.15_dp,'NLME rate')
  call check(sum(nfit%random_effects(:,1)*true_b)>0.0_dp,'NLME random effects direction')

  call acf_values(gfit%residuals,4,acf,pairs,status)
  call check(status==NLME_SUCCESS,'ACF status')
  call check_close(acf(0),1.0_dp,1.0e-12_dp,'ACF lag zero')
  call empirical_variogram(gfit%residuals,xd(:,1),2.0_dp,vg,status)
  call check(status==NLME_SUCCESS .and. sum(vg%pairs)>0,'variogram')

  xs(:,1)=1.0_dp;xs(:,2)=xnm(:,1);zs(:,1)=1.0_dp;gmat(1,1)=0.04_dp
  call simulate_lme(xs,zs,groups,[1.0_dp,0.5_dp],gmat,0.1_dp,ysim,status,seed=123)
  call check(status==NLME_SUCCESS .and. size(ysim)==nn,'simulate LME')
  call check(pooled_sd(ysim,groups,status)>0.0_dp .and. status==NLME_SUCCESS,'pooled SD')
  call check(is_balanced(groups,xnm(:,1)),'balanced grouped data')
  call group_summary(ysim,groups,levels,means,sds,counts,status)
  call check(status==NLME_SUCCESS .and. all(counts==m),'group summary')
  call fit_lm_list(ysim,xs,groups,levels,coefs,sigmas,status)
  call check(status==NLME_SUCCESS .and. size(coefs,1)==ng,'lmList')
  theta0=[1.7_dp,0.5_dp,0.2_dp]
  call fit_nls_list(exponential_decay_model,ynm,xnm,groups,theta0,levels,nlspars,nlssigma,status,control=ctl)
  call check(status==NLME_SUCCESS .and. size(nlspars,1)==ng,'nlsList')
  write(*,'(a)')'test_nonlinear_diagnostics: PASS'
end program test_nonlinear_diagnostics
