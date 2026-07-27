! SPDX-License-Identifier: GPL-3.0-only
! Derived from HDShOP 0.1.7.
module hdshop_inference
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hdshop_kinds, only: dp
  use hdshop_linalg, only: inverse_matrix, quadratic_form
  use hdshop_stats, only: row_means, sample_covariance, normal_quantile, normal_cdf
  use hdshop_portfolio, only: q_matrix, alpha_mv_lt
  implicit none
  private

  type, public :: mvsp_test_result
    real(dp) :: alpha_hat=0.0_dp, alpha_sd=0.0_dp
    real(dp) :: alpha_lower=0.0_dp, alpha_upper=0.0_dp
    real(dp) :: statistic=0.0_dp, p_value=1.0_dp
    logical :: ok=.false.
    character(len=160) :: message=''
  end type mvsp_test_result

  public :: test_mvsp, omega_hat_alpha, d0_vector
  public :: bop19_s, bop19_r, bop19_sigma_s2, bop19_omega

contains

  function d0_vector(gamma,p,n) result(d)
    real(dp),intent(in)::gamma
    integer,intent(in)::p,n
    real(dp)::d(5),ig,c
    c=real(p,dp)/real(n,dp)
    if(ieee_is_finite(gamma))then;ig=1.0_dp/gamma;else;ig=0.0_dp;endif
    d=[(1.0_dp+1.0_dp/(1.0_dp-c))*ig,-1.0_dp,ig*ig/(1.0_dp-c), &
       (-1.0_dp-1.0_dp/(1.0_dp-c))*ig,1.0_dp]
  end function d0_vector

  function omega_hat_alpha(x,b,ok) result(m)
    real(dp),intent(in)::x(:,:),b(:)
    logical,intent(out),optional::ok
    real(dp)::m(5,5)
    real(dp),allocatable::cov(:,:),inv(:,:),mu(:),ones(:),q(:,:)
    real(dp)::c,den,vc,vb,sc,rb,rg,diff
    logical::good
    integer::p,n
    p=size(x,1);n=size(x,2);m=0.0_dp
    cov=sample_covariance(x);call inverse_matrix(cov,inv,good)
    if(.not.good .or. p>=n .or. size(b)/=p)then;if(present(ok))ok=.false.;return;endif
    mu=row_means(x);allocate(ones(p));ones=1.0_dp;c=real(p,dp)/real(n,dp)
    den=dot_product(ones,matmul(inv,ones));vc=(1.0_dp/den)/(1.0_dp-c)
    vb=quadratic_form(b,cov);q=q_matrix(inv);sc=(1.0_dp-c)*quadratic_form(mu,q)-c
    rb=dot_product(b,mu);rg=dot_product(matmul(inv,ones),mu)/den;diff=rb-rg
    m(1,1)=vc*(sc+1.0_dp)/(1.0_dp-c);m(2,2)=2.0_dp*vc*vc/(1.0_dp-c)
    m(3,3)=2.0_dp*((sc+1.0_dp)**2+c-1.0_dp)/(1.0_dp-c)
    m(4,4)=vb;m(5,5)=2.0_dp*vb*vb
    m(4,1)=vc;m(1,4)=vc;m(5,1)=-2.0_dp*vc*diff;m(1,5)=m(5,1)
    m(5,2)=2.0_dp*vc*vc;m(2,5)=m(5,2);m(4,3)=2.0_dp*diff;m(3,4)=m(4,3)
    m(5,3)=-2.0_dp*diff*diff;m(3,5)=m(5,3)
    if(present(ok))ok=.true.
  end function omega_hat_alpha

  function test_mvsp(gamma,x,w0,beta) result(res)
    real(dp),intent(in)::gamma,x(:,:),w0(:)
    real(dp),intent(in),optional::beta
    type(mvsp_test_result)::res
    real(dp)::m(5,5),d(5),level,z,bhat,c,ig,den,vc,vb,sc,rg,rb
    real(dp),allocatable::cov(:,:),inv(:,:),mu(:),ones(:),q(:,:)
    logical::ok
    integer::p,n
    p=size(x,1);n=size(x,2);level=0.05_dp;if(present(beta))level=beta
    if(p>=n .or. size(w0)/=p)then;res%message='test_mvsp requires p < n and matching weights';return;endif
    m=omega_hat_alpha(x,w0,ok);if(.not.ok)then;res%message='covariance inversion failed';return;endif
    cov=sample_covariance(x);call inverse_matrix(cov,inv,ok);mu=row_means(x);allocate(ones(p));ones=1.0_dp
    c=real(p,dp)/real(n,dp);den=dot_product(ones,matmul(inv,ones));vc=(1.0_dp/den)/(1.0_dp-c)
    q=q_matrix(inv);sc=(1.0_dp-c)*quadratic_form(mu,q)-c;rg=dot_product(matmul(inv,ones),mu)/den
    rb=dot_product(w0,mu);vb=quadratic_form(w0,cov)
    if(ieee_is_finite(gamma))then;ig=1.0_dp/gamma;else;ig=0.0_dp;endif
    bhat=vc/(1.0_dp-c)-2.0_dp*(vc+(rb-rg)*ig/(1.0_dp-c))+ &
      (sc+c)*ig*ig/(1.0_dp-c)**3+vb
    res%alpha_hat=alpha_mv_lt(gamma,c,sc,rg,rb,vc,vb)
    d=d0_vector(gamma,p,n)
    res%alpha_sd=sqrt(max(dot_product(d,matmul(m,d)),0.0_dp))/abs(bhat)/sqrt(real(n,dp))
    if(res%alpha_sd<=0.0_dp)then;res%message='zero standard error';return;endif
    z=normal_quantile(1.0_dp-level/2.0_dp);res%alpha_lower=res%alpha_hat-z*res%alpha_sd
    res%alpha_upper=res%alpha_hat+z*res%alpha_sd;res%statistic=res%alpha_hat/res%alpha_sd
    res%p_value=2.0_dp*(1.0_dp-normal_cdf(abs(res%statistic)));res%ok=.true.
  end function test_mvsp

  real(dp) function bop19_s(mu0,inv,mu) result(value)
    real(dp),intent(in)::mu0(:),inv(:,:),mu(:)
    value=quadratic_form(mu,inv)-dot_product(mu0,matmul(inv,mu))**2/quadratic_form(mu0,inv)
  end function bop19_s

  real(dp) function bop19_r(mu0,inv,mu) result(value)
    real(dp),intent(in)::mu0(:),inv(:,:),mu(:)
    value=dot_product(mu0,matmul(inv,mu))/quadratic_form(mu0,inv)
  end function bop19_r

  pure real(dp) function bop19_sigma_s2(c,s) result(value)
    real(dp),intent(in)::c,s
    value=2.0_dp*(c+2.0_dp*s)+2.0_dp*(c+s)**2/(1.0_dp-c)
  end function bop19_sigma_s2

  function bop19_omega(c,sigma_s2,inv,mu0,s,r) result(omega)
    real(dp),intent(in)::c,sigma_s2,inv(:,:),mu0(:),s,r
    real(dp)::omega(2,2),base
    base=c*c*sigma_s2/(c+s)**4
    omega(1,1)=base;omega(1,2)=base*r;omega(2,1)=omega(1,2)
    omega(2,2)=base*r*r+c*c/(c+s)**2*(1.0_dp+(c+s)/(1.0_dp-c))/quadratic_form(mu0,inv)
  end function bop19_omega

end module hdshop_inference
