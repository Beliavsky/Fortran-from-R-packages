! SPDX-License-Identifier: Artistic-2.0
module ecld_models
  use ecd_kinds, only : dp, pi, sqrt_pi, ecd_ok, ecd_invalid
  use ecd_math, only : nan_dp, integrate_adaptive, brent_root, upper_incomplete_gamma, &
    hypergeom_2f0, gamma_p, erfq
  use ecd_core, only : ecd_model, ecd_solve, ecd_new
  implicit none
  private

  type, public :: ecld_model
    real(dp) :: lambda=3.0_dp, sigma=1.0_dp, beta=0.0_dp, mu=0.0_dp
    real(dp) :: epsilon=0.0_dp, rho=0.0_dp, mu_d=0.0_dp
    logical :: is_sged=.false., has_mu_d=.false.
  end type ecld_model

  public :: ecld_new, ecld_from_sd, ecld_solve, ecld_laplace_b
  public :: ecld_const, ecld_pdf, ecld_cdf, ecld_ccdf
  public :: ecld_moment, ecld_mean, ecld_variance, ecld_sd, ecld_skewness, ecld_kurtosis
  public :: ecld_mgf, ecld_imgf, ecld_ogf, ecld_mu_d
  public :: ecld_mgf_quartic, ecld_imgf_quartic, ecld_ogf_quartic
  public :: ecld_mu_d_quartic
  public :: ecld_y_slope, ecld_y_slope_trunc
  public :: ecld_gamma, ecld_gamma_hgeo, ecld_gamma_2f0
  public :: ecld_op_o, ecld_op_v, ecld_op_q, ecld_op_q_skew, ecld_op_u_lag
  public :: ecld_op_vl_quartic, ecld_quartic_q, ecld_quartic_qp
  public :: ecld_quartic_qp_atm_ki, ecld_quartic_qp_rho
  public :: ecld_quartic_qp_skew, ecld_quartic_qp_atm_skew
  public :: ecld_quartic_sn0_atm_ki, ecld_quartic_sn0_rho_stdev
  public :: ecld_quartic_sn0_skew, ecld_quartic_sn0_max_rnv
  public :: ecld_quartic_model_sample
  public :: ecld_fixed_point_atm_ki, ecld_fixed_point_shift
  public :: ecld_incomplete_moment, ecld_imnt_sum
  public :: ecld_ogf_star, ecld_ogf_star_hgeo, ecld_ogf_star_exp
  public :: ecld_ogf_star_gamma_star, ecld_ogf_star_analytic
  public :: ecld_sged_const, ecld_sged_cdf, ecld_sged_moment
  public :: ecld_sged_mgf, ecld_sged_imgf, ecld_sged_ogf

contains

  function ecld_new(lambda,sigma,beta,mu,epsilon,rho,is_sged) result(d)
    real(dp), intent(in), optional :: lambda,sigma,beta,mu,epsilon,rho
    logical, intent(in), optional :: is_sged
    type(ecld_model) :: d
    if(present(lambda)) d%lambda=lambda
    if(present(sigma)) d%sigma=sigma
    if(present(beta)) d%beta=beta
    if(present(mu)) d%mu=mu
    if(present(epsilon)) d%epsilon=epsilon
    if(present(rho)) d%rho=rho
    if(present(is_sged)) d%is_sged=is_sged
  end function ecld_new

  function ecld_from_sd(lambda,sd,beta,mu,status) result(d)
    real(dp), intent(in), optional :: lambda,sd,beta,mu
    integer, intent(out), optional :: status
    type(ecld_model) :: d
    real(dp) :: lam,sdev,b,s0,lo,hi
    integer :: st
    lam=3.0_dp; sdev=1.0_dp; b=0.0_dp
    if(present(lambda))lam=lambda
    if(present(sd))sdev=sd
    if(present(beta))b=beta
    d=ecld_new(lam,1.0_dp,b)
    if(present(mu))d%mu=mu
    if(present(status))status=ecd_ok
    if(sdev<=0.0_dp .or. lam<=0.0_dp) then
      if(present(status))status=ecd_invalid; d%sigma=nan_dp(); return
    end if
    s0=sdev*sqrt(gamma(lam/2.0_dp)/gamma(1.5_dp*lam))
    if(b==0.0_dp) then
      d%sigma=s0
    else
      lo=0.1_dp*s0; hi=10.0_dp*s0
      d%sigma=brent_root(match_sd,lo,hi,1e-10_dp,100,st)
      if(present(status))status=st
    end if
    d%mu_d=ecld_mu_d(d); d%has_mu_d=.true.
  contains
    function match_sd(s) result(v)
      real(dp), intent(in) :: s
      real(dp) :: v
      type(ecld_model) :: q
      q=ecld_new(lam,s,b)
      v=ecld_sd(q)-sdev
    end function match_sd
  end function ecld_from_sd

  function ecld_solve(d,x) result(y)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: y,xi,s2,b0
    type(ecd_model) :: e
    xi=(x-d%mu)/d%sigma
    if(d%is_sged) then
      s2=d%sigma*merge(1.0_dp-d%beta,1.0_dp+d%beta,xi<0.0_dp)
      y=-abs((x-d%mu)/s2)**(2.0_dp/d%lambda)
    else if(d%beta==0.0_dp) then
      y=-abs(xi)**(2.0_dp/d%lambda)
    else if(d%lambda==2.0_dp) then
      b0=ecld_laplace_b(d%beta,merge(1.0_dp,-1.0_dp,xi<0.0_dp))
      y=-b0*abs(xi)
    else
      e%alpha=0.0_dp; e%gamma=0.0_dp; e%sigma=d%sigma; e%beta=d%beta
      e%mu=d%mu; e%lambda=d%lambda; e%cusp=0; e%has_const=.false.
      y=ecd_solve(e,x)
    end if
  end function ecld_solve

  pure elemental function ecld_laplace_b(beta,sgn,sigma) result(v)
    real(dp), intent(in) :: beta
    real(dp), intent(in), optional :: sgn,sigma
    real(dp) :: v,g,s
    g=0.0_dp; s=0.0_dp
    if(present(sgn))g=sgn
    if(present(sigma))s=sigma
    v=sqrt(1.0_dp+beta*beta/4.0_dp)+g*beta/2.0_dp+g*s
  end function ecld_laplace_b

  function ecld_const(d) result(c)
    type(ecld_model), intent(in) :: d
    real(dp) :: c
    if(d%is_sged .or. d%beta==0.0_dp) then
      c=d%sigma*d%lambda*gamma(d%lambda/2.0_dp)
    else if(d%lambda==2.0_dp) then
      c=2.0_dp*d%sigma*sqrt(1.0_dp+d%beta*d%beta/4.0_dp)
    else
      c=integrate_adaptive(core,-huge(1.0_dp),huge(1.0_dp),1e-9_dp,1e-11_dp)
    end if
  contains
    function core(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v,y
      y=ecld_solve(d,x)
      if(y<log(tiny(1.0_dp))) then; v=0.0_dp; else; v=exp(y); end if
    end function core
  end function ecld_const

  function ecld_pdf(d,x) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v,y
    y=ecld_solve(d,x)
    if(y<log(tiny(1.0_dp))) then; v=0.0_dp; else; v=exp(y)/ecld_const(d); end if
  end function ecld_pdf

  function ecld_cdf(d,x) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v
    if(.not.d%is_sged .and. d%beta==0.0_dp) then
      if(x==d%mu) then
        v=0.5_dp
      else if(x<d%mu) then
        v=0.5_dp*(1.0_dp-gamma_p(d%lambda/2.0_dp,abs((x-d%mu)/d%sigma)**(2.0_dp/d%lambda)))
      else
        v=0.5_dp+0.5_dp*gamma_p(d%lambda/2.0_dp,abs((x-d%mu)/d%sigma)**(2.0_dp/d%lambda))
      end if
      return
    end if
    if(x<d%mu) then
      v=integrate_adaptive(pdf_fun,-huge(1.0_dp),x,1e-9_dp,1e-11_dp)
    else
      v=1.0_dp-integrate_adaptive(pdf_fun,x,huge(1.0_dp),1e-9_dp,1e-11_dp)
    end if
    v=max(0.0_dp,min(1.0_dp,v))
  contains
    function pdf_fun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=ecld_pdf(d,z)
    end function pdf_fun
  end function ecld_cdf

  function ecld_ccdf(d,x) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v
    v=1.0_dp-ecld_cdf(d,x)
  end function ecld_ccdf

  function ecld_moment(d,order) result(v)
    type(ecld_model), intent(in) :: d
    integer, intent(in) :: order
    real(dp) :: v,sn,sp,b0,bp,bn
    if(d%is_sged) then
      sn=d%sigma*(1.0_dp-d%beta); sp=d%sigma*(1.0_dp+d%beta)
      v=(sp**(order+1)+(-1.0_dp)**order*sn**(order+1))* &
        gamma(d%lambda*real(order+1,dp)/2.0_dp)/(2.0_dp*d%sigma*gamma(d%lambda/2.0_dp))
    else if(d%beta==0.0_dp) then
      if(mod(order,2)/=0) then; v=0.0_dp
      else; v=d%sigma**order*gamma(d%lambda*real(order+1,dp)/2.0_dp)/gamma(d%lambda/2.0_dp); end if
    else if(d%lambda==2.0_dp) then
      b0=ecld_laplace_b(d%beta); bp=ecld_laplace_b(d%beta,1.0_dp); bn=ecld_laplace_b(d%beta,-1.0_dp)
      v=d%sigma**order*gamma(real(order+1,dp))/(2.0_dp*b0)* &
        (bp**(order+1)+(-1.0_dp)**order*bn**(order+1))
    else
      v=integrate_adaptive(mfun,-huge(1.0_dp),huge(1.0_dp),1e-8_dp,1e-10_dp)
    end if
  contains
    function mfun(x) result(q)
      real(dp), intent(in) :: x
      real(dp) :: q
      q=(x-d%mu)**order*ecld_pdf(d,x)
    end function mfun
  end function ecld_moment

  function ecld_mean(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v
    v=d%mu+ecld_moment(d,1)
  end function ecld_mean

  function ecld_variance(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v,m1
    m1=ecld_moment(d,1); v=ecld_moment(d,2)-m1*m1
  end function ecld_variance

  function ecld_sd(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v
    v=sqrt(max(0.0_dp,ecld_variance(d)))
  end function ecld_sd

  function ecld_skewness(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v,m1,m2,m3,var
    m1=ecld_moment(d,1); m2=ecld_moment(d,2); m3=ecld_moment(d,3); var=m2-m1*m1
    v=(m3-3*m1*m2+2*m1**3)/var**1.5_dp
  end function ecld_skewness

  function ecld_kurtosis(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v,m1,m2,m3,m4,var
    if(.not.d%is_sged .and. d%beta==0.0_dp) then
      v=gamma(d%lambda/2.0_dp)*gamma(2.5_dp*d%lambda)/gamma(1.5_dp*d%lambda)**2
      return
    end if
    m1=ecld_moment(d,1); m2=ecld_moment(d,2); m3=ecld_moment(d,3); m4=ecld_moment(d,4)
    var=m2-m1*m1
    v=(m4-4*m1*m3+6*m1*m1*m2-3*m1**4)/var**2
  end function ecld_kurtosis

  function ecld_mgf(d,t) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in), optional :: t
    real(dp) :: v,tt,term,last
    integer :: n
    tt=1.0_dp; if(present(t))tt=t
    if(.not.d%is_sged .and. d%lambda==4.0_dp .and. d%beta==0.0_dp) then
      v=ecld_mgf_quartic(d,tt); return
    end if
    if(.not.d%is_sged .and. d%lambda==1.0_dp .and. d%beta==0.0_dp) then
      v=exp(tt*d%mu+d%sigma*d%sigma*tt*tt/4.0_dp); return
    end if
    if(.not.d%is_sged .and. d%lambda==2.0_dp) then
      v=exp(tt*d%mu)/(1.0_dp-d%beta*d%sigma*tt-d%sigma*d%sigma*tt*tt); return
    end if
    if(d%is_sged .and. d%lambda==2.0_dp) then
      v=exp(tt*d%mu)/(1.0_dp-2*d%beta*d%sigma*tt-(1-d%beta*d%beta)*d%sigma*d%sigma*tt*tt); return
    end if
    v=1.0_dp; last=huge(1.0_dp)
    do n=1,500
      term=ecld_moment(d,n)*tt**n/gamma(real(n+1,dp))
      if(abs(term)>last .and. n>10) exit
      v=v+term; last=abs(term)
      if(abs(term)<1e-13_dp*max(1.0_dp,abs(v))) exit
    end do
    v=exp(tt*d%mu)*v
  end function ecld_mgf

  function ecld_imgf(d,k,option_type,risk_neutral) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    logical, intent(in), optional :: risk_neutral
    real(dp) :: v,shift
    character(len=1) :: ot
    logical :: rn
    ot='c'; if(present(option_type))ot=option_type(1:1)
    rn=.true.; if(present(risk_neutral))rn=risk_neutral
    shift=0.0_dp
    if(rn) then
      if(.not.d%has_mu_d) then; shift=ecld_mu_d(d)-d%mu; else; shift=d%mu_d-d%mu; end if
    end if
    if(.not.d%is_sged .and. d%lambda==4.0_dp .and. d%beta==0.0_dp) then
      v=ecld_imgf_quartic(d,k,ot,rn); return
    end if
    if(ot=='c' .or. ot=='C') then
      v=integrate_adaptive(ifun,k,huge(1.0_dp),1e-8_dp,1e-11_dp)
    else
      v=integrate_adaptive(ifun,-huge(1.0_dp),k,1e-8_dp,1e-11_dp)
    end if
  contains
    function ifun(x) result(q)
      real(dp), intent(in) :: x
      real(dp) :: q
      q=exp(x)*ecld_pdf_shifted(d,x,shift)
    end function ifun
  end function ecld_imgf

  function ecld_ogf(d,k,option_type,risk_neutral) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    logical, intent(in), optional :: risk_neutral
    real(dp) :: v,shift
    character(len=1) :: ot
    logical :: rn
    ot='c'; if(present(option_type))ot=option_type(1:1)
    rn=.true.; if(present(risk_neutral))rn=risk_neutral
    shift=0.0_dp
    if(rn) then
      if(.not.d%has_mu_d) then; shift=ecld_mu_d(d)-d%mu; else; shift=d%mu_d-d%mu; end if
    end if
    if(.not.d%is_sged .and. d%lambda==4.0_dp .and. d%beta==0.0_dp) then
      v=ecld_ogf_quartic(d,k,ot,rn); return
    end if
    if(ot=='c' .or. ot=='C') then
      v=integrate_adaptive(cfun,k,huge(1.0_dp),1e-8_dp,1e-11_dp)
    else
      v=integrate_adaptive(pfun,-huge(1.0_dp),k,1e-8_dp,1e-11_dp)
    end if
  contains
    function cfun(x) result(q)
      real(dp), intent(in) :: x
      real(dp) :: q
      q=(exp(x)-exp(k))*ecld_pdf_shifted(d,x,shift)
    end function cfun
    function pfun(x) result(q)
      real(dp), intent(in) :: x
      real(dp) :: q
      q=(exp(k)-exp(x))*ecld_pdf_shifted(d,x,shift)
    end function pfun
  end function ecld_ogf

  function ecld_mgf_quartic(d,t) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in), optional :: t
    real(dp) :: v,tt,st,z,dz,mz
    tt=1.0_dp;if(present(t))tt=t
    if(d%lambda/=4.0_dp .or. d%beta/=0.0_dp .or. d%sigma*tt<=0.0_dp)then
      v=nan_dp();return
    end if
    st=d%sigma*tt
    z=1.0_dp/sqrt(4.0_dp*st)
    dz=z*z*exp(-z*z)
    mz=z**3*(erfq(z,-1)-erfq(z,1))
    v=exp(tt*d%mu)*(mz+dz)
  end function ecld_mgf_quartic

  function ecld_mu_d_quartic(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v,z,dz,mz
    if(d%lambda/=4.0_dp .or. d%beta/=0.0_dp .or. d%sigma<=0.0_dp)then
      v=nan_dp();return
    end if
    z=1.0_dp/(2.0_dp*sqrt(d%sigma))
    dz=z*z*exp(-z*z)
    mz=z**3*(erfq(z,-1)-erfq(z,1))
    v=-log(mz+dz)
  end function ecld_mu_d_quartic

  function ecld_imgf_quartic(d,k,option_type,risk_neutral) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    logical, intent(in), optional :: risk_neutral
    real(dp) :: v,s,mu,m1,ki,ki_inf,z,dz,mc_plus,mp_minus
    character(len=1) :: ot
    logical :: rn
    s=d%sigma;ot='c';rn=.true.
    if(present(option_type))ot=option_type(1:1)
    if(present(risk_neutral))rn=risk_neutral
    if(d%lambda/=4.0_dp .or. d%beta/=0.0_dp .or. s<=0.0_dp)then
      v=nan_dp();return
    end if
    if(rn)then
      if(d%has_mu_d)then;mu=d%mu_d;else;mu=ecld_mu_d_quartic(d);end if
      m1=1.0_dp
    else
      mu=d%mu;m1=ecld_mgf_quartic(d)
    end if
    ki=(k-mu)/s
    ki_inf=1.0_dp/(4.0_dp*s*s)
    ki=min(ki,ki_inf)
    z=1.0_dp/(2.0_dp*sqrt(s));dz=z*z*exp(-z*z)
    mc_plus=exp(mu)*(quartic_ms(-1.0_dp,ki,s,z)+dz)
    mp_minus=exp(mu)*quartic_ms(1.0_dp,ki,s,z)
    if(ot=='p'.or.ot=='P')then
      if(ki<0.0_dp)then;v=mp_minus;else;v=m1-mc_plus;end if
    else
      if(ki>=0.0_dp)then;v=mc_plus;else;v=m1-mp_minus;end if
    end if
  end function ecld_imgf_quartic

  function quartic_ms(sgn,ki,s,z) result(v)
    real(dp), intent(in) :: sgn,ki,s,z
    real(dp) :: v,abki,vv,w
    abki=abs(ki)
    vv=sqrt(abki)+sgn*abki*s
    w=z+sgn*sqrt(abki*s)
    v=exp(-vv)*(z**3*erfq(w,nint(sgn))-z*z)
  end function quartic_ms

  function ecld_ogf_quartic(d,k,option_type,risk_neutral) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    logical, intent(in), optional :: risk_neutral
    real(dp) :: v,s,mu,ki,ki2,z,dz,mz,m0kz,lp_raw,lc_raw
    character(len=1) :: ot
    logical :: rn
    s=d%sigma;ot='c';rn=.true.
    if(present(option_type))ot=option_type(1:1)
    if(present(risk_neutral))rn=risk_neutral
    if(d%lambda/=4.0_dp .or. d%beta/=0.0_dp .or. s<=0.0_dp)then
      v=nan_dp();return
    end if
    if(rn)then
      if(d%has_mu_d)then;mu=d%mu_d;else;mu=ecld_mu_d_quartic(d);end if
    else
      mu=d%mu
    end if
    ki=(k-mu)/s;ki2=sqrt(abs(ki));z=1.0_dp/(2.0_dp*sqrt(s))
    dz=z*z*exp(-z*z);mz=z**3*(erfq(z,-1)-erfq(z,1))
    m0kz=mz-exp(s*ki)+dz
    lp_raw=exp(-ki2*(1.0_dp+s*ki2))* &
      (-z*z+z**3*erfq(z+ki2/(2.0_dp*z),1)+0.5_dp*(ki2+1.0_dp))
    lc_raw=exp(-ki2*(1.0_dp-s*ki2))* &
      (-z*z+z**3*erfq(z-ki2/(2.0_dp*z),-1)-0.5_dp*(ki2+1.0_dp))+dz
    if(ot=='p'.or.ot=='P')then
      if(ki<0.0_dp)then;v=exp(mu)*lp_raw;else;v=exp(mu)*(lc_raw-m0kz);end if
    else
      if(ki>=0.0_dp)then;v=exp(mu)*lc_raw;else;v=exp(mu)*(lp_raw+m0kz);end if
    end if
  end function ecld_ogf_quartic

  function ecld_pdf_shifted(d,x,shift) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: x,shift
    real(dp) :: v
    type(ecld_model) :: q
    q=d; q%mu=d%mu+shift
    v=ecld_pdf(q,x)
  end function ecld_pdf_shifted

  function ecld_mu_d(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v
    v=-log(ecld_mgf(d,1.0_dp))
  end function ecld_mu_d

  function ecld_y_slope(d,x) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v,xi,y,s2
    xi=(x-d%mu)/d%sigma; y=ecld_solve(d,x)
    if(d%is_sged) then
      s2=d%sigma*merge(1.0_dp-d%beta,1.0_dp+d%beta,xi<0.0_dp)
      if(x==d%mu) then; v=0.0_dp; else
        v=-(2.0_dp/d%lambda)*sign(abs((x-d%mu)/s2)**(2.0_dp/d%lambda-1.0_dp),x-d%mu)/s2
      end if
    else
      v=-(d%beta*y+2.0_dp*xi)/(d%lambda*(-y)**(d%lambda-1.0_dp)+d%beta*xi)/d%sigma
    end if
  end function ecld_y_slope

  function ecld_y_slope_trunc(d,t) result(x)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in), optional :: t
    real(dp) :: x,tt,lo,hi
    integer :: st
    tt=1.0_dp; if(present(t))tt=t
    lo=d%mu; hi=d%mu+max(10.0_dp*d%sigma,1.0_dp)
    do while(ecld_y_slope(d,hi)+tt>0.0_dp .and. hi<d%mu+1e6_dp*d%sigma); hi=d%mu+2*(hi-d%mu); end do
    x=brent_root(sf,lo+epsilon(1.0_dp),hi,1e-10_dp,200,st)
  contains
    function sf(z) result(v)
      real(dp), intent(in) :: z
      real(dp) :: v
      v=ecld_y_slope(d,z)+tt
    end function sf
  end function ecld_y_slope_trunc

  function ecld_gamma(s,x) result(v)
    real(dp), intent(in) :: s
    real(dp), intent(in), optional :: x
    real(dp) :: v,z
    z=0.0_dp; if(present(x))z=x
    v=upper_incomplete_gamma(s,z)
  end function ecld_gamma

  function ecld_gamma_hgeo(s,x,order) result(v)
    real(dp), intent(in) :: s,x
    integer, intent(in) :: order
    real(dp) :: v
    v=exp(-x)*x**(s-1.0_dp)*hypergeom_2f0(s,x,order)
  end function ecld_gamma_hgeo

  function ecld_gamma_2f0(s,x,order) result(v)
    real(dp), intent(in) :: s,x
    integer, intent(in) :: order
    real(dp) :: v
    v=hypergeom_2f0(s,x,order)
  end function ecld_gamma_2f0

  function ecld_op_o(sigma1,k,option_type,rho) result(l)
    real(dp), intent(in) :: sigma1,k
    character(len=*), intent(in), optional :: option_type
    real(dp), intent(in), optional :: rho
    real(dp) :: l,r,p,q,m1k,sgn
    character(len=1) :: ot
    r=0.0_dp; if(present(rho))r=rho
    ot='c'; if(present(option_type))ot=option_type(1:1)
    if(sigma1<=0.0_dp) then; l=nan_dp(); return; end if
    m1k=exp(r)-exp(k)
    p=-0.5_dp*exp(r)*erf((k-r)/sigma1-sigma1/4.0_dp)
    q=0.5_dp*exp(k)*erf((k-r)/sigma1+sigma1/4.0_dp)
    sgn=merge(1.0_dp,-1.0_dp,ot=='c' .or. ot=='C')
    l=p+q+sgn*m1k/2.0_dp
  end function ecld_op_o

  function ecld_op_v(l,k,option_type,rho,ttm,status) result(sigma1)
    real(dp), intent(in) :: l,k
    character(len=*), intent(in), optional :: option_type
    real(dp), intent(in), optional :: rho,ttm
    integer, intent(out), optional :: status
    real(dp) :: sigma1,r,lo,hi,tm
    character(len=1) :: ot
    integer :: st
    r=0.0_dp; if(present(rho))r=rho
    ot='c'; if(present(option_type))ot=option_type(1:1)
    if(present(status))status=ecd_ok
    if(l<=0.0_dp) then; sigma1=nan_dp(); if(present(status))status=ecd_invalid; return; end if
    lo=1e-8_dp; hi=max(0.1_dp,100.0_dp*l)
    do while(rootf(hi)<0.0_dp .and. hi<100.0_dp); hi=hi*2.0_dp; end do
    sigma1=brent_root(rootf,lo,hi,1e-10_dp,300,st)
    if(present(status))status=st
    if(present(ttm)) then
      tm=ttm; if(tm>0.0_dp)sigma1=sigma1/sqrt(2.0_dp*tm)
    end if
  contains
    function rootf(s) result(v)
      real(dp), intent(in) :: s
      real(dp) :: v
      v=ecld_op_o(s,k,ot,r)-l
    end function rootf
  end function ecld_op_v

  function ecld_op_q(d,ki,option_type) result(q)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    character(len=*), intent(in), optional :: option_type
    real(dp) :: q,k,l
    character(len=1) :: ot
    ot='c'; if(present(option_type))ot=option_type(1:1)
    k=ki*d%sigma+d%mu
    l=ecld_ogf(d,k,ot,.false.)+d%epsilon
    q=ecld_op_v(l,k,ot)/d%sigma
  end function ecld_op_q

  function ecld_op_q_skew(d,ki,dki,option_type) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    real(dp), intent(in), optional :: dki
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v,h
    h=0.1_dp; if(present(dki))h=dki
    v=(ecld_op_q(d,ki+h/2.0_dp,option_type)-ecld_op_q(d,ki-h/2.0_dp,option_type))/h
  end function ecld_op_q_skew

  subroutine ecld_op_u_lag(local_prices,log_strikes,scale,n_lag,slope,status)
    real(dp), intent(in) :: local_prices(:),log_strikes(:),scale
    integer, intent(in), optional :: n_lag
    real(dp), allocatable, intent(out) :: slope(:)
    integer, intent(out), optional :: status
    integer :: n,i,m
    if(present(status))status=ecd_ok
    n=2;if(present(n_lag))n=n_lag
    if(size(local_prices)/=size(log_strikes) .or. n<1 .or. n>=size(local_prices) .or. any(local_prices<=0.0_dp))then
      allocate(slope(0));if(present(status))status=ecd_invalid;return
    end if
    m=size(local_prices)-n
    allocate(slope(m))
    do i=1,m
      slope(i)=(log(local_prices(i+n))-log(local_prices(i)))/ &
        (log_strikes(i+n)-log_strikes(i))*scale
    end do
  end subroutine ecld_op_u_lag

  function ecld_op_vl_quartic(d,k,option_type,ttm,status) result(iv)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    real(dp), intent(in), optional :: ttm
    integer, intent(out), optional :: status
    real(dp) :: iv,l,ks
    character(len=1) :: ot
    integer :: st
    ot='c';if(present(option_type))ot=option_type(1:1)
    if(abs(d%lambda-4.0_dp)>1.0e-12_dp .or. d%beta/=0.0_dp)then
      iv=nan_dp();if(present(status))status=ecd_invalid;return
    end if
    ks=k-d%rho
    l=ecld_ogf_quartic(d,ks,ot,.false.)+d%epsilon
    if(present(ttm))then
      iv=ecld_op_v(l,ks,ot,0.0_dp,ttm,st)
    else
      iv=ecld_op_v(l,ks,ot,0.0_dp,status=st)
    end if
    if(present(status))status=st
  end function ecld_op_vl_quartic

  function ecld_quartic_q(d,ki,option_type,status) result(q)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    character(len=*), intent(in), optional :: option_type
    integer, intent(out), optional :: status
    real(dp) :: q,k,iv
    integer :: st
    k=ki*d%sigma+d%mu+d%rho
    iv=ecld_op_vl_quartic(d,k,option_type,status=st)
    q=iv/d%sigma
    if(present(status))status=st
  end function ecld_quartic_q

  function ecld_quartic_qp(d,ki,status) result(q)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(out), optional :: status
    real(dp) :: q
    q=ecld_quartic_q(d,ki,'p',status)
  end function ecld_quartic_qp

  function ecld_quartic_qp_atm_ki(d,lower,upper,status) result(root)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in), optional :: lower,upper
    integer, intent(out), optional :: status
    real(dp) :: root,lo,hi,flo,fhi
    integer :: st,n
    lo=-50.0_dp;hi=-1.37_dp
    if(present(lower))lo=lower
    if(present(upper))hi=upper
    flo=target(lo);fhi=target(hi);n=0
    do while(flo*fhi>=0.0_dp .and. n<60)
      if(flo<0.0_dp)then;hi=lo;fhi=flo;lo=1.5_dp*lo;flo=target(lo)
      else;hi=hi/1.5_dp;fhi=target(hi);end if
      n=n+1
    end do
    if(flo*fhi>=0.0_dp)then;root=nan_dp();st=ecd_invalid
    else;root=brent_root(target,lo,hi,1.0e-9_dp,200,st);end if
    if(present(status))status=st
  contains
    function target(ki) result(v)
      real(dp), intent(in)::ki
      real(dp)::v
      v=ecld_quartic_qp(d,ki)-sqrt(240.0_dp)
    end function target
  end function ecld_quartic_qp_atm_ki

  function ecld_quartic_qp_rho(d,atm_ki,lower,upper,status) result(rho)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in), optional :: atm_ki,lower,upper
    integer, intent(out), optional :: status
    real(dp) :: rho,a
    integer :: st
    if(present(atm_ki))then;a=atm_ki;st=ecd_ok
    else;a=ecld_quartic_qp_atm_ki(d,lower,upper,st);end if
    rho=-(d%mu+a*d%sigma)
    if(present(status))status=st
  end function ecld_quartic_qp_rho

  function ecld_quartic_qp_skew(d,ki,dki) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    real(dp), intent(in), optional :: dki
    real(dp) :: v,h
    h=0.1_dp;if(present(dki))h=dki
    v=(ecld_quartic_qp(d,ki+h)-ecld_quartic_qp(d,ki))/h
  end function ecld_quartic_qp_skew

  function ecld_quartic_qp_atm_skew(d,dki,lower,upper,status) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in), optional :: dki,lower,upper
    integer, intent(out), optional :: status
    real(dp) :: v,a,h
    integer :: st
    h=0.1_dp;if(present(dki))h=dki
    a=ecld_quartic_qp_atm_ki(d,lower,upper,st)
    v=ecld_quartic_qp_skew(d,a,h)
    if(present(status))status=st
  end function ecld_quartic_qp_atm_skew

  function ecld_quartic_sn0_atm_ki(status) result(root)
    integer, intent(out), optional :: status
    real(dp) :: root,q0
    integer :: st
    q0=sqrt(240.0_dp)
    root=brent_root(eq,-30.0_dp,0.0_dp,1.0e-11_dp,200,st)
    if(present(status))status=st
  contains
    function eq(ki) result(v)
      real(dp),intent(in)::ki
      real(dp)::v
      v=q0*ecld_ogf_star_analytic(ecld_new(lambda=1.0_dp),abs(ki/q0))- &
        ecld_ogf_star_analytic(ecld_new(lambda=4.0_dp),abs(ki))
    end function eq
  end function ecld_quartic_sn0_atm_ki

  function ecld_quartic_sn0_rho_stdev() result(v)
    real(dp)::v
    v=-ecld_quartic_sn0_atm_ki()/sqrt(120.0_dp)
  end function ecld_quartic_sn0_rho_stdev

  function ecld_quartic_sn0_skew() result(v)
    real(dp)::v,ki,x,z
    ki=ecld_quartic_sn0_atm_ki();x=sqrt(abs(ki));z=abs(ki/sqrt(240.0_dp))
    v=sqrt(pi)*exp(z*z-x)*(1.0_dp+x-exp(x)*erfc(z))
  end function ecld_quartic_sn0_skew

  function quartic_skew_root(mu_ratio,status) result(atm)
    real(dp), intent(in) :: mu_ratio
    integer, intent(out), optional :: status
    real(dp) :: atm
    integer :: st
    atm=brent_root(skew_eq,-4.0_dp,-0.1_dp,1.0e-8_dp,100,st)
    if(present(status))status=st
  contains
    function skew_eq(ki) result(zv)
      real(dp),intent(in)::ki
      real(dp)::zv,x,z
      x=sqrt(abs(ki))
      z=abs(ki/sqrt(240.0_dp)-mu_ratio*sqrt(120.0_dp)/sqrt(240.0_dp))
      zv=sqrt(pi)*exp(z*z-x)*(1.0_dp+x-exp(x)*erfc(z))
    end function skew_eq
  end function quartic_skew_root

  function quartic_q_root(mu_ratio,atm,status) result(q)
    real(dp), intent(in) :: mu_ratio,atm
    integer, intent(out), optional :: status
    real(dp) :: q
    integer :: st
    q=brent_root(q_eq,0.1_dp,100.0_dp,1.0e-8_dp,100,st)
    if(present(status))status=st
  contains
    function q_eq(qv) result(zv)
      real(dp),intent(in)::qv
      real(dp)::zv,ki1
      ki1=atm/qv-mu_ratio*sqrt(120.0_dp)/qv
      zv=qv*ecld_ogf_star_analytic(ecld_new(lambda=1.0_dp),ki1)- &
        ecld_ogf_star_analytic(ecld_new(lambda=4.0_dp),atm)
    end function q_eq
  end function quartic_q_root

  function ecld_quartic_sn0_max_rnv(status) result(root)
    integer, intent(out), optional :: status
    real(dp)::root
    integer::st
    root=brent_root(outer,0.29_dp,0.30_dp,1.0e-9_dp,100,st)
    if(present(status))status=st
  contains
    function outer(mu_ratio) result(v)
      real(dp),intent(in)::mu_ratio
      real(dp)::v,atm,q
      integer::st1,st2
      atm=quartic_skew_root(mu_ratio,st1)
      if(st1/=ecd_ok)then
        v=huge(1.0_dp)
        return
      end if
      q=quartic_q_root(mu_ratio,atm,st2)
      if(st2/=ecd_ok)then
        v=huge(1.0_dp)
        return
      end if
      v=(q-sqrt(240.0_dp))*1000.0_dp
    end function outer
  end function ecld_quartic_sn0_max_rnv

  subroutine ecld_quartic_model_sample(date_code,ttm,sigma,mu_plus_ratio,epsilon_ratio,skew_adjusted,status)
    character(len=*),intent(in)::date_code
    real(dp),intent(in)::ttm(:)
    real(dp),intent(out)::sigma(:),mu_plus_ratio(:),epsilon_ratio(:)
    logical,intent(in),optional::skew_adjusted
    integer,intent(out),optional::status
    logical::sa
    sa=.true.;if(present(skew_adjusted))sa=skew_adjusted
    if(size(sigma)/=size(ttm).or.size(mu_plus_ratio)/=size(ttm).or.size(epsilon_ratio)/=size(ttm))then
      if(present(status))status=ecd_invalid;return
    end if
    select case(trim(date_code))
    case('2015-07-20')
      sigma=0.08_dp*(1.0_dp+sqrt(ttm))*sqrt(ttm/120.0_dp);epsilon_ratio=0.0136_dp
      if(sa)then;mu_plus_ratio=-(0.230_dp-0.53_dp*sqrt(ttm))
      else;mu_plus_ratio=-(0.29_dp-0.479_dp*sqrt(ttm));end if
    case('2015-08-24')
      sigma=0.2_dp*ttm**(-0.27_dp)*sqrt(ttm/120.0_dp);epsilon_ratio=0.0174_dp
      if(sa)then;mu_plus_ratio=-(0.315_dp-0.21_dp*sqrt(ttm))
      else;mu_plus_ratio=-(0.295_dp-0.185_dp*sqrt(ttm));end if
    case('2015-10-06')
      sigma=0.18_dp*sqrt(ttm/120.0_dp);epsilon_ratio=0.0094_dp
      if(sa)then;mu_plus_ratio=-(0.325_dp-0.44_dp*sqrt(ttm))
      else;mu_plus_ratio=-(0.293_dp-0.337_dp*sqrt(ttm));end if
    case default
      sigma=nan_dp();mu_plus_ratio=nan_dp();epsilon_ratio=nan_dp()
      if(present(status))status=ecd_invalid;return
    end select
    if(present(status))status=ecd_ok
  end subroutine ecld_quartic_model_sample

  function ecld_incomplete_moment(d,ki,order,option_type) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(in) :: order
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v,k,m1,g2,k2,scale_side,kii,nsgn
    character(len=1) :: ot
    ot='c'; if(present(option_type))ot=option_type(1:1)
    if(order<0) then; v=nan_dp(); return; end if
    if((.not.d%is_sged .and. d%beta==0.0_dp) .or. d%is_sged) then
      scale_side=d%sigma
      kii=ki
      if(d%is_sged) then
        scale_side=d%sigma*merge(1.0_dp-d%beta,1.0_dp+d%beta,ki<0.0_dp)
        kii=ki*d%sigma/scale_side
      end if
      k2=abs(kii)**(2.0_dp/d%lambda)
      nsgn=merge((-1.0_dp)**order,1.0_dp,ki<0.0_dp)
      g2=nsgn*upper_incomplete_gamma(d%lambda*real(order+1,dp)/2.0_dp,k2)/ &
        gamma(d%lambda/2.0_dp)/2.0_dp*scale_side**order
      if(d%is_sged)g2=g2*scale_side/d%sigma
      m1=ecld_moment(d,order)
      if(ot=='c'.or.ot=='C') then; v=merge(g2,m1-g2,ki>=0.0_dp)
      else; v=merge(g2,m1-g2,ki<0.0_dp); end if
      return
    end if
    k=d%mu+d%sigma*ki
    if(ot=='c'.or.ot=='C') then
      v=integrate_adaptive(integrand,k,huge(1.0_dp),1e-8_dp,1e-11_dp)
    else
      v=integrate_adaptive(integrand,-huge(1.0_dp),k,1e-8_dp,1e-11_dp)
    end if
  contains
    function integrand(x) result(q)
      real(dp), intent(in) :: x
      real(dp) :: q
      q=(x-d%mu)**order*ecld_pdf(d,x)
    end function integrand
  end function ecld_incomplete_moment

  function ecld_imnt_sum(d,ki,max_order,option_type) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(in) :: max_order
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v
    integer :: n
    v=0.0_dp
    do n=0,max(0,min(max_order,100))
      v=v+ecld_incomplete_moment(d,ki,n,option_type)/gamma(real(n+1,dp))
    end do
  end function ecld_imnt_sum

  function ecld_ogf_star(d,ki) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    real(dp) :: v,xi
    xi=abs(ki)**(2.0_dp/d%lambda)
    v=(upper_incomplete_gamma(d%lambda,xi)-abs(ki)* &
      upper_incomplete_gamma(d%lambda/2.0_dp,xi))/(2.0_dp*gamma(d%lambda/2.0_dp))
  end function ecld_ogf_star

  function ecld_ogf_star_hgeo(d,ki,order) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(in), optional :: order
    real(dp) :: v,xi,q
    integer :: n
    n=4; if(present(order))n=order
    if(ki==0.0_dp)then;v=ecld_ogf_star(d,ki);return;end if
    xi=abs(ki)**(2.0_dp/d%lambda)
    q=exp(-xi)*abs(ki)**(2.0_dp-2.0_dp/d%lambda)/(2.0_dp*gamma(d%lambda/2.0_dp))
    v=q*(hypergeom_2f0(d%lambda,xi,n)-hypergeom_2f0(d%lambda/2.0_dp,xi,n))
  end function ecld_ogf_star_hgeo

  function ecld_ogf_star_exp(d,ki,order) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(in), optional :: order
    real(dp) :: v,xi,p,l,dfac
    integer :: n,i
    n=3; if(present(order))n=order
    if(ki==0.0_dp)then;v=ecld_ogf_star(d,ki);return;end if
    l=d%lambda; xi=abs(ki)**(2.0_dp/l)
    v=exp(-xi)*abs(ki)**(2.0_dp-4.0_dp/l)*l/(4.0_dp*gamma(l/2.0_dp))
    if(n==0)return
    p=1.0_dp
    if(l==1.0_dp)then
      dfac=1.0_dp
      do i=1,n;dfac=dfac*real(2*i+1,dp);p=p+dfac/(-2.0_dp*xi)**i;end do
      v=v*p;return
    end if
    if(n>=1)p=p+1.5_dp*(l-2.0_dp)/xi
    if(n>=2)p=p+0.25_dp*(7*l*l-36*l+44)/xi**2
    if(n>=3)p=p+0.125_dp*(15*l**3-140*l*l+420*l-400)/xi**3
    v=v*p
  end function ecld_ogf_star_exp

  function ecld_ogf_star_gamma_star(d,ki,order) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(in), optional :: order
    real(dp) :: v,k2l,q,p
    integer :: n,nmax
    nmax=6;if(present(order))nmax=order
    k2l=abs(ki)**(2.0_dp/d%lambda)
    q=gamma(d%lambda)/(2.0_dp*gamma(d%lambda/2.0_dp))-abs(ki)/2.0_dp
    p=0.0_dp
    do n=0,nmax
      p=p+k2l**n*(1.0_dp/gamma(d%lambda/2.0_dp+real(n+1,dp))- &
        gamma(d%lambda)/(gamma(d%lambda/2.0_dp)*gamma(d%lambda+real(n+1,dp))))
    end do
    v=q+0.5_dp*ki*ki*exp(-k2l)*p
  end function ecld_ogf_star_gamma_star

  function ecld_ogf_star_analytic(d,ki,status) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: ki
    integer, intent(out), optional :: status
    real(dp) :: v,x
    x=abs(ki);if(present(status))status=ecd_ok
    select case(nint(d%lambda))
    case(1)
      if(abs(d%lambda-1.0_dp)>1e-12_dp)then;v=nan_dp();if(present(status))status=ecd_invalid;return;end if
      v=exp(-x*x)/(2.0_dp*sqrt(pi))-0.5_dp*x*erfc(x)
    case(2)
      if(abs(d%lambda-2.0_dp)>1e-12_dp)then;v=nan_dp();if(present(status))status=ecd_invalid;return;end if
      v=0.5_dp*exp(-x)
    case(3)
      if(abs(d%lambda-3.0_dp)>1e-12_dp)then;v=nan_dp();if(present(status))status=ecd_invalid;return;end if
      v=2.0_dp/sqrt(pi)*exp(-x**(2.0_dp/3.0_dp))*(1.0_dp+x**(2.0_dp/3.0_dp))- &
        0.5_dp*x*erfc(x**(1.0_dp/3.0_dp))
    case(4)
      if(abs(d%lambda-4.0_dp)>1e-12_dp)then;v=nan_dp();if(present(status))status=ecd_invalid;return;end if
      v=exp(-sqrt(x))*(3.0_dp+3.0_dp*sqrt(x)+x)
    case default
      v=nan_dp();if(present(status))status=ecd_invalid
    end select
  end function ecld_ogf_star_analytic

  function ecld_sged_const(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v
    v=ecld_const(d)
  end function ecld_sged_const
  function ecld_sged_cdf(d,x) result(v)
    type(ecld_model), intent(in) :: d; real(dp), intent(in) :: x; real(dp) :: v
    v=ecld_cdf(d,x)
  end function ecld_sged_cdf
  function ecld_sged_moment(d,n) result(v)
    type(ecld_model), intent(in) :: d; integer, intent(in) :: n; real(dp) :: v
    v=ecld_moment(d,n)
  end function ecld_sged_moment
  function ecld_sged_mgf(d,t) result(v)
    type(ecld_model), intent(in) :: d; real(dp), intent(in), optional :: t; real(dp) :: v
    if(present(t))then;v=ecld_mgf(d,t);else;v=ecld_mgf(d);end if
  end function ecld_sged_mgf
  function ecld_sged_imgf(d,k,t,option_type) result(v)
    type(ecld_model), intent(in) :: d; real(dp), intent(in) :: k
    real(dp), intent(in), optional :: t
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v,tt
    type(ecld_model) :: q
    tt=1.0_dp;if(present(t))tt=t
    q=d;q%sigma=d%sigma*tt;q%mu=d%mu*tt
    v=ecld_imgf(q,tt*k,option_type,.false.)
  end function ecld_sged_imgf
  function ecld_sged_ogf(d,k,option_type) result(v)
    type(ecld_model), intent(in) :: d; real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v
    v=ecld_ogf(d,k,option_type,.false.)
  end function ecld_sged_ogf

  pure function ecld_fixed_point_atm_ki(d) result(v)
    type(ecld_model), intent(in) :: d
    real(dp) :: v
    v=(-d%mu-d%rho)/d%sigma
  end function ecld_fixed_point_atm_ki

  pure function ecld_fixed_point_shift(d,atm_imp_k) result(v)
    type(ecld_model), intent(in) :: d
    real(dp), intent(in) :: atm_imp_k
    real(dp) :: v
    v=-(atm_imp_k-d%mu)
  end function ecld_fixed_point_shift


end module ecld_models
