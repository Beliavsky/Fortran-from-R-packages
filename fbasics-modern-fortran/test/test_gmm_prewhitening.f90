! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module gmm_prewhite_test_moments
  use fbasics_kinds, only: dp
  implicit none
contains
  subroutine serial_regression_moments(theta,data,g)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::g(:,:)
    real(dp),allocatable::e(:)
    integer::n
    n=size(data,1);allocate(e(n),g(n,3))
    e=data(:,2)-theta(1)-theta(2)*data(:,1)
    g(:,1)=e
    g(:,2)=e*data(:,1)
    g(:,3)=e*(data(:,1)**2-sum(data(:,1)**2)/real(n,dp))
  end subroutine serial_regression_moments
end module gmm_prewhite_test_moments

program test_gmm_prewhitening
  use fbasics
  use test_support
  use gmm_prewhite_test_moments
  implicit none
  integer,parameter::n_arma=520,n_var=640,n_reg=260
  integer::t,info
  real(dp),parameter::rho_true=0.58_dp,psi_true=0.42_dp
  real(dp)::rho,psi,sigma,obj,bw_arma,bw_ar1,bw_expected,alpha2,common,denom,bw_used
  real(dp),allocatable::arma(:),innov(:),gm(:,:),g(:,:),eps(:,:),resid(:,:),recolor(:,:),coef(:,:,:)
  real(dp),allocatable::gc(:,:),manual(:,:),s_pw(:,:),s_direct(:,:),system(:,:),expected_recolor(:,:)
  real(dp),allocatable::data(:,:)
  real(dp)::a(2,2),determinant
  logical::ok,pw_ok
  type(gmm_result)::fit

  call set_lcg_seed(202607249_8)

  allocate(arma(n_arma),innov(n_arma),gm(n_arma,1))
  do t=1,n_arma
    innov(t)=rnorm_fs(0.0_dp,1.0_dp)
  end do
  arma(1)=innov(1)
  do t=2,n_arma
    arma(t)=rho_true*arma(t-1)+innov(t)+psi_true*innov(t-1)
  end do
  call fit_arma11_css(arma,rho,psi,sigma,ok,obj,max_iter=420)
  call assert_true(ok,'ARMA(1,1) CSS convergence')
  call assert_close(rho,rho_true,0.10_dp,'ARMA(1,1) AR coefficient')
  call assert_close(psi,psi_true,0.12_dp,'ARMA(1,1) MA coefficient')
  call assert_true(sigma>0.75_dp.and.sigma<1.25_dp,'ARMA(1,1) innovation scale')
  gm(:,1)=arma
  bw_arma=andrews_bandwidth_value(gm,'Quadratic Spectral','ARMA(1,1)',ok)
  call assert_true(ok.and.bw_arma>=1.0_dp,'ARMA(1,1) Andrews bandwidth')
  call fit_arma11_css(arma-sum(arma)/real(n_arma,dp),rho,psi,sigma,ok,max_iter=420)
  denom=((1.0_dp+psi)*sigma/(1.0_dp-rho))**4
  common=((1.0_dp+rho*psi)*(rho+psi))**2*sigma**4
  alpha2=4.0_dp*common/(1.0_dp-rho)**8/denom
  bw_expected=1.3221_dp*(real(n_arma,dp)*alpha2)**0.2_dp
  call assert_close(bw_arma,bw_expected,2.0e-7_dp,'ARMA(1,1) Andrews formula')
  bw_ar1=andrews_bandwidth_value(gm,'Quadratic Spectral','AR(1)',ok)
  call assert_true(ok.and.abs(bw_arma-bw_ar1)>0.05_dp,'ARMA and AR bandwidth paths differ')

  allocate(g(n_var,2),eps(n_var,2))
  a=reshape([0.55_dp,-0.08_dp,0.10_dp,0.35_dp],[2,2])
  do t=1,n_var
    eps(t,1)=rnorm_fs(0.0_dp,1.0_dp)
    eps(t,2)=rnorm_fs(0.0_dp,0.8_dp)
  end do
  g(1,:)=eps(1,:)
  do t=2,n_var
    g(t,:)=matmul(a,g(t-1,:))+eps(t,:)
  end do
  call prewhiten_var(g,1,resid,recolor,pw_ok,coef)
  call assert_true(pw_ok.and.size(resid,1)==n_var-1,'VAR(1) prewhitening')
  call assert_close(coef(1,1,1),a(1,1),0.08_dp,'VAR prewhite A11')
  call assert_close(coef(1,2,1),a(1,2),0.08_dp,'VAR prewhite A12')
  call assert_close(coef(2,1,1),a(2,1),0.08_dp,'VAR prewhite A21')
  call assert_close(coef(2,2,1),a(2,2),0.08_dp,'VAR prewhite A22')
  allocate(system(2,2));system=0.0_dp;system(1,1)=1.0_dp;system(2,2)=1.0_dp
  system=system-coef(:,:,1)
  call matrix_inverse(system,expected_recolor,info)
  call assert_true(info==0,'VAR recoloring inverse')
  call assert_close(maxval(abs(recolor-expected_recolor)),0.0_dp,2.0e-11_dp,'VAR recoloring matrix')
  call prewhiten_var(g,2,resid,recolor,pw_ok,coef)
  call assert_true(pw_ok.and.size(coef,3)==2,'higher-order VAR prewhitening')
  call assert_true(maxval(abs(coef(:,:,2)))<0.13_dp,'redundant second VAR lag')

  allocate(gc(n_var,2));gc=g-spread(sum(g,dim=1)/real(n_var,dp),1,n_var)
  call prewhiten_var(gc,1,resid,recolor,pw_ok)
  manual=matmul(recolor,matmul(matmul(transpose(resid),resid)/real(n_var,dp),transpose(recolor)))
  call moment_covariance(g,s_pw,cov_type='HAC',kernel='Bartlett',bandwidth=0,prewhite_order=1, &
    bandwidth_method='fixed',bandwidth_used=bw_used,prewhite_succeeded=pw_ok)
  call assert_true(pw_ok.and.abs(bw_used)<1.0e-14_dp,'prewhitened fixed-bandwidth HAC')
  call assert_close(maxval(abs(s_pw-manual)),0.0_dp,2.0e-10_dp,'prewhitening recoloring identity')
  call moment_covariance(g,s_direct,cov_type='HAC',kernel='Bartlett',bandwidth=0,prewhite_order=0,bandwidth_method='fixed')
  call assert_true(maxval(abs(s_pw-s_direct))>0.05_dp,'prewhitening changes serial covariance estimate')
  determinant=s_pw(1,1)*s_pw(2,2)-s_pw(1,2)*s_pw(2,1)
  call assert_true(all(s_pw==s_pw).and.minval([s_pw(1,1),s_pw(2,2),determinant])>-1.0e-10_dp,'prewhitened HAC positive semidefinite')

  allocate(data(n_reg,3));innov=0.0_dp
  do t=1,n_reg
    data(t,1)=-1.0_dp+2.0_dp*real(t-1,dp)/real(n_reg-1,dp)
    if(t==1)then
      innov(t)=rnorm_fs(0.0_dp,1.0_dp)
      data(t,3)=innov(t)
    else
      innov(t)=rnorm_fs(0.0_dp,1.0_dp)
      data(t,3)=0.55_dp*data(t-1,3)+innov(t)+0.30_dp*innov(t-1)
    end if
    data(t,2)=1.0_dp+2.0_dp*data(t,1)+0.08_dp*data(t,3)
  end do
  call fit_gmm(serial_regression_moments,data,[0.0_dp,0.0_dp],[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp], &
    'two_step',fit,cov_type='HAC',kernel='Quadratic Spectral',bandwidth=-1,max_iter=500, &
    prewhite_order=1,bandwidth_method='andrews_arma11')
  call assert_close(fit%theta(1),1.0_dp,0.06_dp,'prewhitened GMM intercept')
  call assert_close(fit%theta(2),2.0_dp,0.08_dp,'prewhitened GMM slope')
  call assert_true(fit%prewhite_succeeded.and.fit%prewhite_order==1,'GMM prewhitening status')
  call assert_true(trim(fit%bandwidth_method)=='andrews_arma11'.and.fit%bandwidth_used>=1.0_dp,'GMM ARMA bandwidth status')
  call assert_true(all(fit%covariance==fit%covariance).and.maxval(abs(fit%covariance-transpose(fit%covariance)))<1.0e-9_dp,'prewhitened GMM covariance')

  g=1.0_dp
  call moment_covariance(g,s_pw,cov_type='HAC',kernel='Quadratic Spectral',bandwidth=-1,prewhite_order=1, &
    bandwidth_method='andrews_arma11',bandwidth_used=bw_used,prewhite_succeeded=pw_ok)
  call assert_true(.not.pw_ok.and.all(s_pw==s_pw),'singular prewhitening fallback')

  write(*,'(a)')'GMM prewhitening and Andrews ARMA(1,1) bandwidth tests passed.'
end program test_gmm_prewhitening
