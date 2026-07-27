! SPDX-License-Identifier: GPL-3.0-only
! Scalar and matrix formulas translated from HDShOP 0.1.7.
module hdshop_formulas
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use hdshop_kinds, only: dp
  use hdshop_linalg, only: inverse_matrix, quadratic_form
  use hdshop_portfolio, only: q_matrix, r_gmv, v_gmv, v_portfolio
  implicit none
  private
  public :: alpha_star, alpha_star_gmv, var_alpha_simple
  public :: alpha_star_n_bop19, beta_star_n_bop19
  public :: alpha_star_bop19, beta_star_bop19
  public :: omega_lest, omega_lest_old

contains

  pure real(dp) function reciprocal_gamma(gamma) result(value)
    real(dp),intent(in)::gamma
    if(ieee_is_finite(gamma))then;value=1.0_dp/gamma;else;value=0.0_dp;endif
  end function reciprocal_gamma

  real(dp) function alpha_star(gamma,mu,covariance,b,c) result(alpha)
    real(dp),intent(in)::gamma,mu(:),covariance(:,:),b(:),c
    real(dp),allocatable::inv(:,:),q(:,:)
    real(dp)::rg,rb,vg,vb,ss,ig,num,den
    logical::ok
    call inverse_matrix(covariance,inv,ok)
    if(.not.ok)then;alpha=huge(1.0_dp);return;endif
    q=q_matrix(inv);rg=r_gmv(mu,covariance);rb=dot_product(b,mu)
    vg=v_gmv(covariance);vb=v_portfolio(covariance,b);ss=quadratic_form(mu,q);ig=reciprocal_gamma(gamma)
    num=(rg-rb)*(1.0_dp+1.0_dp/(1.0_dp-c))*ig+(vb-vg)+ss*ig*ig/(1.0_dp-c)
    den=vg/(1.0_dp-c)-2.0_dp*(vg+(rb-rg)*ig/(1.0_dp-c))+ &
      (ss+c)*ig*ig/(1.0_dp-c)**3+vb
    alpha=num/den
  end function alpha_star

  real(dp) function alpha_star_gmv(covariance,b,c) result(alpha)
    real(dp),intent(in)::covariance(:,:),b(:),c
    real(dp)::vg,vb,num
    vg=v_gmv(covariance);vb=v_portfolio(covariance,b);num=(1.0_dp-c)*(vb-vg)
    alpha=num/(num+c*vg)
  end function alpha_star_gmv

  real(dp) function var_alpha_simple(covariance,b,n) result(value)
    real(dp),intent(in)::covariance(:,:),b(:)
    integer,intent(in)::n
    real(dp)::c,vb,vg,lb
    c=real(size(covariance,1),dp)/real(n,dp);vb=v_portfolio(covariance,b);vg=v_gmv(covariance)
    lb=vb/vg-1.0_dp
    value=2.0_dp*(1.0_dp-c)*c*c*(lb+1.0_dp)*((2.0_dp-c)*lb+c)/ &
      ((1.0_dp-c)*lb+c)**4
  end function var_alpha_simple

  real(dp) function alpha_star_n_bop19(y,inv,mu,mu0) result(alpha)
    real(dp),intent(in)::y(:),inv(:,:),mu(:),mu0(:)
    real(dp)::y_mu,y_0,mu_0,zero_0,y_y,den
    y_mu=dot_product(y,matmul(inv,mu));y_0=dot_product(y,matmul(inv,mu0))
    mu_0=dot_product(mu,matmul(inv,mu0));zero_0=quadratic_form(mu0,inv);y_y=quadratic_form(y,inv)
    den=y_y*zero_0-y_0*y_0;alpha=(y_mu*zero_0-mu_0*y_0)/den
  end function alpha_star_n_bop19

  real(dp) function beta_star_n_bop19(y,inv,mu,mu0) result(beta)
    real(dp),intent(in)::y(:),inv(:,:),mu(:),mu0(:)
    real(dp)::y_y,y_0,y_mu,mu_0,zero_0,den
    y_y=quadratic_form(y,inv);y_0=dot_product(y,matmul(inv,mu0));y_mu=dot_product(y,matmul(inv,mu))
    mu_0=dot_product(mu,matmul(inv,mu0));zero_0=quadratic_form(mu0,inv)
    den=y_y*zero_0-y_0*y_0;beta=(y_y*mu_0-y_0*y_mu)/den
  end function beta_star_n_bop19

  real(dp) function alpha_star_bop19(c,mu,inv,mu0) result(alpha)
    real(dp),intent(in)::c,mu(:),inv(:,:),mu0(:)
    real(dp)::i1,i2,i3
    i1=quadratic_form(mu,inv);i2=quadratic_form(mu0,inv);i3=dot_product(mu,matmul(inv,mu0))
    alpha=(i1*i2-i3*i3)/((c+i1)*i2-i3*i3)
  end function alpha_star_bop19

  real(dp) function beta_star_bop19(alpha,mu,inv,mu0) result(beta)
    real(dp),intent(in)::alpha,mu(:),inv(:,:),mu0(:)
    beta=(1.0_dp-alpha)*dot_product(mu,matmul(inv,mu0))/quadratic_form(mu0,inv)
  end function beta_star_bop19

  function omega_lest(s,c,gamma,v_c,l,q,eta) result(omega)
    real(dp),intent(in)::s,c,gamma,v_c,l(:,:),q(:,:),eta(:)
    real(dp),allocatable::omega(:,:)
    real(dp)::ig,a,b
    integer::m
    m=size(l,1);allocate(omega(m,m));ig=reciprocal_gamma(gamma)
    a=(ig*ig*(s+1.0_dp)+v_c)*(1.0_dp-c)
    b=ig*ig*(s+c)**2
    omega=a*matmul(matmul(l,q),transpose(l))+b*spread(eta,2,m)*spread(eta,1,m)
  end function omega_lest

  function omega_lest_old(s,c,gamma,v_c,l,q,eta) result(omega)
    real(dp),intent(in)::s,c,gamma,v_c,l(:,:),q(:,:),eta(:)
    real(dp),allocatable::omega(:,:)
    real(dp)::ig,a,b
    integer::m
    m=size(l,1);allocate(omega(m,m));ig=reciprocal_gamma(gamma)
    a=(((1.0_dp-c)/(s+c)+(s+c)*ig)*ig+v_c)*(1.0_dp-c)
    b=ig*ig*(2.0_dp*(1.0_dp-c)*c**3/(s+c)**2+ &
      4.0_dp*(1.0_dp-c)*c*s*(s+2.0_dp*c)/(s+c)**2+ &
      2.0_dp*(1.0_dp-c)*c*c*(s+c)**2/(s*s)-s*s)
    omega=a*matmul(matmul(l,q),transpose(l))+b*spread(eta,2,m)*spread(eta,1,m)
  end function omega_lest_old

end module hdshop_formulas
