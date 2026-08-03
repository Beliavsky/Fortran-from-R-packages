! SPDX-License-Identifier: Artistic-2.0
module ecd_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ecd_kinds, only : dp, pi, sqrt_pi, ecd_ok, ecd_invalid
  use ecd_math, only : nan_dp, cbrt_real, integrate_adaptive, brent_root, gamma_p, gamma_quantile
  use ecd_rng, only : rng_state, rng_uniform
  implicit none
  private

  type, public :: ecd_stats_type
    real(dp) :: m1=0.0_dp, m2=0.0_dp, m3=0.0_dp, m4=0.0_dp
    real(dp) :: mean=0.0_dp, variance=0.0_dp, stdev=0.0_dp
    real(dp) :: skewness=0.0_dp, kurtosis=0.0_dp
  end type ecd_stats_type

  type, public :: ecd_model
    real(dp) :: alpha=0.0_dp, gamma=0.0_dp, sigma=1.0_dp
    real(dp) :: beta=0.0_dp, mu=0.0_dp, lambda=3.0_dp
    integer :: cusp=1
    real(dp) :: radius=0.0_dp, theta=0.0_dp
    real(dp) :: norm_const=1.5_dp*sqrt_pi
    logical :: has_const=.true., has_stats=.false.
    type(ecd_stats_type) :: stats
  end type ecd_model

  type, public :: ellipticity_result
    real(dp) :: left=0.0_dp, right=0.0_dp, average=0.0_dp
  end type ellipticity_result

  public :: ecd_new, ecd_cusp_new, ecd_polar_new, ecd_initialize
  public :: ecd_cusp_a2r, ecd_cusp_r2a, ecd_adj_gamma, ecd_adj2gamma
  public :: ecd_discriminant, ecd_j_invariant, ecd_y0_isomorphic
  public :: ecd_solve, ecd_solve_sym, ecd_solve_trig, ecd_solve_cusp_asym
  public :: ecd_y_slope, ecd_pdf, ecd_cdf, ecd_ccdf, ecd_quantile, ecd_random
  public :: ecd_moment, ecd_statistics, ecd_asymptotic_statistics, ecd_ellipticity
  public :: ecd_cusp_std_moment, ecd_cusp_std_cf, ecd_cusp_std_mgf
  public :: ecd_imgf, ecd_mu_d, ecd_ogf, ecd_rational, ecd_max_kurtosis

contains

  function ecd_new(alpha,gamma,sigma,beta,mu,lambda,with_stats,status) result(d)
    real(dp), intent(in), optional :: alpha,gamma,sigma,beta,mu,lambda
    logical, intent(in), optional :: with_stats
    integer, intent(out), optional :: status
    type(ecd_model) :: d
    logical :: ws
    if(present(alpha)) d%alpha=alpha
    if(present(gamma)) d%gamma=gamma
    if(present(sigma)) d%sigma=sigma
    if(present(beta)) d%beta=beta
    if(present(mu)) d%mu=mu
    if(present(lambda)) d%lambda=lambda
    d%cusp=0
    if(d%alpha==0.0_dp .and. d%gamma==0.0_dp .and. d%beta==0.0_dp .and. d%lambda==3.0_dp) d%cusp=1
    ws=.false.; if(present(with_stats)) ws=with_stats
    call ecd_initialize(d,ws,status)
  end function ecd_new

  function ecd_cusp_new(alpha,gamma,sigma,mu,with_stats,status) result(d)
    real(dp), intent(in), optional :: alpha,gamma,sigma,mu
    logical, intent(in), optional :: with_stats
    integer, intent(out), optional :: status
    type(ecd_model) :: d
    logical :: ws
    d%sigma=1.0_dp; d%mu=0.0_dp
    if(present(sigma)) d%sigma=sigma
    if(present(mu)) d%mu=mu
    if(present(alpha)) then
      d%alpha=alpha; d%gamma=ecd_cusp_a2r(alpha); d%cusp=1
    else if(present(gamma)) then
      d%gamma=gamma; d%alpha=ecd_cusp_r2a(gamma); d%cusp=2
    else
      d%alpha=0.0_dp; d%gamma=0.0_dp; d%cusp=1
    end if
    ws=.false.; if(present(with_stats)) ws=with_stats
    call ecd_initialize(d,ws,status)
  end function ecd_cusp_new

  function ecd_polar_new(radius,theta,sigma,beta,mu,with_stats,status) result(d)
    real(dp), intent(in) :: radius,theta
    real(dp), intent(in), optional :: sigma,beta,mu
    logical, intent(in), optional :: with_stats
    integer, intent(out), optional :: status
    type(ecd_model) :: d
    real(dp) :: ag
    logical :: ws
    d%radius=radius; d%theta=theta
    d%alpha=radius*cos(theta)
    ag=radius*sin(theta)
    d%gamma=ecd_adj2gamma(ag)
    if(present(sigma)) d%sigma=sigma
    if(present(beta)) d%beta=beta
    if(present(mu)) d%mu=mu
    ws=.false.; if(present(with_stats)) ws=with_stats
    call ecd_initialize(d,ws,status)
  end function ecd_polar_new

  subroutine ecd_initialize(d,with_stats,status)
    type(ecd_model), intent(inout) :: d
    logical, intent(in), optional :: with_stats
    integer, intent(out), optional :: status
    real(dp) :: ag
    logical :: ws
    if(present(status)) status=ecd_ok
    if(d%sigma<=0.0_dp .or. d%lambda<=0.0_dp) then
      if(present(status)) status=ecd_invalid
      d%norm_const=nan_dp(); d%has_const=.false.; return
    end if
    ag=ecd_adj_gamma(d%gamma)
    d%radius=sqrt(d%alpha*d%alpha+ag*ag)
    if(d%radius==0.0_dp) then
      d%theta=0.0_dp
    else
      d%theta=acos(max(-1.0_dp,min(1.0_dp,d%alpha/d%radius)))
      if(ag<0.0_dp) d%theta=2.0_dp*pi-d%theta
    end if
    if(d%cusp>0) then
      if(d%alpha<0.0_dp .or. d%gamma>0.0_dp) then
        if(present(status)) status=ecd_invalid
      end if
    end if
    d%norm_const=ecd_normalizing_constant(d)
    d%has_const=ieee_is_finite(d%norm_const) .and. d%norm_const>0.0_dp
    ws=.false.; if(present(with_stats)) ws=with_stats
    if(ws .and. d%has_const) then
      d%stats=ecd_statistics(d)
      d%has_stats=.true.
    end if
  end subroutine ecd_initialize

  pure elemental function ecd_cusp_a2r(alpha) result(gamma)
    real(dp), intent(in) :: alpha
    real(dp) :: gamma
    gamma=-(27.0_dp*abs(alpha)**2/4.0_dp)**(1.0_dp/3.0_dp)
  end function ecd_cusp_a2r

  pure elemental function ecd_cusp_r2a(gamma) result(alpha)
    real(dp), intent(in) :: gamma
    real(dp) :: alpha
    alpha=sqrt(4.0_dp*abs(gamma)**3/27.0_dp)
  end function ecd_cusp_r2a

  pure elemental function ecd_adj_gamma(gamma) result(v)
    real(dp), intent(in) :: gamma
    real(dp) :: v
    v=sign(ecd_cusp_r2a(gamma),gamma)
  end function ecd_adj_gamma

  pure elemental function ecd_adj2gamma(v) result(gamma)
    real(dp), intent(in) :: v
    real(dp) :: gamma
    gamma=sign(abs(ecd_cusp_a2r(v)),v)
  end function ecd_adj2gamma

  pure function ecd_discriminant(d) result(v)
    type(ecd_model), intent(in) :: d
    real(dp) :: v
    v=-16.0_dp*(4.0_dp*d%gamma**3+27.0_dp*d%alpha**2)
  end function ecd_discriminant

  pure function ecd_j_invariant(d) result(v)
    type(ecd_model), intent(in) :: d
    real(dp) :: v,disc
    disc=ecd_discriminant(d)
    if(disc==0.0_dp) then; v=nan_dp(); else; v=-1728.0_dp*(4.0_dp*d%gamma)**3/disc; end if
  end function ecd_j_invariant

  function ecd_solve(d,x) result(y)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: y,xi
    if(d%lambda==3.0_dp) then
      y=ecd_solve_trig(d,x)
      return
    end if
    xi=(x-d%mu)/d%sigma
    if(d%beta==0.0_dp .and. d%gamma==0.0_dp) then
      y=sign(abs(d%alpha-xi*xi)**(1.0_dp/d%lambda),d%alpha-xi*xi)
      return
    end if
    y=general_lambda_root(d,xi)
  end function ecd_solve

  function ecd_solve_sym(d,x) result(y)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: y
    if(d%beta/=0.0_dp) then; y=nan_dp(); else; y=ecd_solve_trig(d,x); end if
  end function ecd_solve_sym

  function ecd_solve_trig(d,x) result(y)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: y,xi,a,g,disc,v,ang,sg,rootarg
    xi=(x-d%mu)/d%sigma
    if(d%lambda/=3.0_dp) then
      y=general_lambda_root(d,xi); return
    end if
    a=d%alpha-xi*xi
    g=d%gamma+d%beta*xi
    if(g==0.0_dp) then
      y=cbrt_real(a); return
    end if
    disc=-16.0_dp*(4.0_dp*g**3+27.0_dp*a*a)
    if(d%cusp>0 .and. abs(disc)<1000.0_dp*epsilon(1.0_dp)) disc=0.0_dp
    sg=sqrt(abs(g)/3.0_dp)
    if(disc>=0.0_dp) then
      v=1.5_dp*a*sqrt(3.0_dp/abs(g))/g
      v=max(-1.0_dp,min(1.0_dp,v))
      ang=acos(v)
      y=-2.0_dp*sg*cos(ang/3.0_dp)
    else if(g<0.0_dp) then
      rootarg=-1.5_dp*abs(a)*sqrt(3.0_dp/abs(g))/g
      rootarg=max(1.0_dp,rootarg)
      y=2.0_dp*sign(1.0_dp,a)*sg*cosh(acosh(rootarg)/3.0_dp)
    else
      v=-1.5_dp*a*sqrt(3.0_dp/abs(g))/g
      y=-2.0_dp*sg*sinh(asinh(v)/3.0_dp)
    end if
  end function ecd_solve_trig

  recursive pure function ecd_solve_cusp_asym(x,beta) result(y)
    real(dp), intent(in) :: x,beta
    real(dp) :: y,x0,v,w,a
    if(beta==0.0_dp) then; y=-abs(x)**(2.0_dp/3.0_dp); return; end if
    if(beta<0.0_dp) then; y=ecd_solve_cusp_asym(-x,-beta); return; end if
    x0=-4.0_dp*beta**3/27.0_dp
    v=sqrt(abs(x/x0)); w=2.0_dp*sqrt(abs(beta*x/3.0_dp))
    if(x>=0.0_dp) then
      y=-w*sinh(asinh(v)/3.0_dp)
    else if(x<x0) then
      y=-w*cosh(acosh(v)/3.0_dp)
    else
      a=acos(max(-1.0_dp,min(1.0_dp,v)))
      y=-w*cos(a/3.0_dp)
    end if
  end function ecd_solve_cusp_asym

  function general_lambda_root(d,xi) result(y)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: xi
    real(dp) :: y,lo,hi,lastx,lastf,x,fx,best
    integer :: i,status
    hi=0.0_dp
    lo=-max(2.0_dp,2.0_dp*abs(xi)**(2.0_dp/d%lambda)+2.0_dp*abs(d%beta)+2.0_dp*abs(d%gamma))
    do while(rootfun(lo)<0.0_dp)
      lo=lo*2.0_dp
      if(abs(lo)>1e12_dp) exit
    end do
    best=nan_dp(); lastx=lo; lastf=rootfun(lastx)
    do i=1,800
      x=lo+(hi-lo)*real(i,dp)/800.0_dp; fx=rootfun(x)
      if(lastf==0.0_dp) best=lastx
      if(lastf*fx<=0.0_dp) then
        best=brent_root(rootfun,lastx,x,1e-12_dp,200,status)
        exit
      end if
      lastx=x; lastf=fx
    end do
    y=best
  contains
    function rootfun(z) result(v)
      real(dp), intent(in) :: z
      real(dp) :: v
      v=(-z)**d%lambda-(d%gamma+d%beta*xi)*z+d%alpha-xi*xi
    end function rootfun
  end function general_lambda_root

  function ecd_y_slope(d,x) result(v)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v,xi,y,den
    xi=(x-d%mu)/d%sigma; y=ecd_solve(d,x)
    den=d%lambda*(-y)**(d%lambda-1.0_dp)+d%beta*xi+d%gamma
    if(den==0.0_dp) then; v=sign(huge(1.0_dp),-(d%beta*y+2.0_dp*xi)); else
      v=-(d%beta*y+2.0_dp*xi)/den/d%sigma
    end if
  end function ecd_y_slope

  function ecd_normalizing_constant(d) result(c)
    type(ecd_model), intent(in) :: d
    real(dp) :: c
    if(d%alpha==0.0_dp .and. d%gamma==0.0_dp .and. d%beta==0.0_dp) then
      c=d%sigma*d%lambda*gamma(d%lambda/2.0_dp)
      return
    end if
    if(d%beta==0.0_dp) then
      c=2.0_dp*d%sigma*integrate_adaptive(core,0.0_dp,huge(1.0_dp),1e-9_dp,1e-11_dp)
    else
      c=d%sigma*integrate_adaptive(core,-huge(1.0_dp),huge(1.0_dp),1e-9_dp,1e-11_dp)
    end if
  contains
    function core(z) result(v)
      real(dp), intent(in) :: z
      real(dp) :: v,y
      y=ecd_solve(d,d%mu+d%sigma*z)
      if(y<log(tiny(1.0_dp))) then; v=0.0_dp; else; v=exp(y); end if
    end function core
  end function ecd_normalizing_constant

  function ecd_pdf(d,x) result(v)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v,y,c
    c=d%norm_const; if(.not.d%has_const) c=ecd_normalizing_constant(d)
    y=ecd_solve(d,x)
    if(y<log(tiny(1.0_dp))) then; v=0.0_dp; else; v=exp(y)/c; end if
  end function ecd_pdf

  function ecd_cdf(d,x,from_x) result(v)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: from_x
    real(dp) :: v,a
    a=-huge(1.0_dp); if(present(from_x)) a=from_x
    if(.not.present(from_x) .and. d%alpha==0.0_dp .and. d%gamma==0.0_dp .and. d%beta==0.0_dp) then
      if(x==d%mu) then
        v=0.5_dp
      else if(x<d%mu) then
        v=0.5_dp*(1.0_dp-gamma_p(d%lambda/2.0_dp,abs((x-d%mu)/d%sigma)**(2.0_dp/d%lambda)))
      else
        v=0.5_dp+0.5_dp*gamma_p(d%lambda/2.0_dp,abs((x-d%mu)/d%sigma)**(2.0_dp/d%lambda))
      end if
      return
    end if
    if(x<=a) then; v=0.0_dp; else; v=integrate_adaptive(pdf_fun,a,x,1e-9_dp,1e-11_dp); end if
    if(.not.present(from_x)) v=max(0.0_dp,min(1.0_dp,v))
  contains
    function pdf_fun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=ecd_pdf(d,z)
    end function pdf_fun
  end function ecd_cdf

  function ecd_ccdf(d,x,to_x) result(v)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: to_x
    real(dp) :: v,b
    b=huge(1.0_dp); if(present(to_x)) b=to_x
    if(x>=b) then; v=0.0_dp; else; v=integrate_adaptive(pdf_fun,x,b,1e-9_dp,1e-11_dp); end if
    if(.not.present(to_x)) v=max(0.0_dp,min(1.0_dp,v))
  contains
    function pdf_fun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=ecd_pdf(d,z)
    end function pdf_fun
  end function ecd_ccdf

  function ecd_quantile(d,p,status) result(x)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: p
    integer, intent(out), optional :: status
    real(dp) :: x,lo,hi,sd
    integer :: st
    if(present(status)) status=ecd_ok
    if(p<0.0_dp .or. p>1.0_dp) then; x=nan_dp(); if(present(status))status=ecd_invalid; return; end if
    if(p==0.0_dp) then; x=-huge(1.0_dp); return; end if
    if(p==1.0_dp) then; x=huge(1.0_dp); return; end if
    if(d%alpha==0.0_dp .and. d%gamma==0.0_dp .and. d%beta==0.0_dp) then
      if(p==0.5_dp) then
        x=d%mu
      else if(p<0.5_dp) then
        x=d%mu-d%sigma*gamma_quantile(1.0_dp-2.0_dp*p,d%lambda/2.0_dp,1.0_dp,st)**(d%lambda/2.0_dp)
      else
        x=d%mu+d%sigma*gamma_quantile(2.0_dp*p-1.0_dp,d%lambda/2.0_dp,1.0_dp,st)**(d%lambda/2.0_dp)
      end if
      if(present(status))status=st
      return
    end if
    if(d%has_stats) then; sd=d%stats%stdev; else; sd=sqrt(max(ecd_moment(d,2,.true.),d%sigma*d%sigma)); end if
    lo=d%mu-8.0_dp*sd; hi=d%mu+8.0_dp*sd
    do while(ecd_cdf(d,lo)>p); lo=d%mu+2.0_dp*(lo-d%mu); end do
    do while(ecd_cdf(d,hi)<p); hi=d%mu+2.0_dp*(hi-d%mu); end do
    x=brent_root(qfun,lo,hi,1e-9_dp,150,st)
    if(present(status)) status=st
  contains
    function qfun(z) result(v)
      real(dp), intent(in) :: z
      real(dp) :: v
      v=ecd_cdf(d,z)-p
    end function qfun
  end function ecd_quantile

  subroutine ecd_random(d,rng,x)
    type(ecd_model), intent(in) :: d
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    integer :: i
    do i=1,size(x); x(i)=ecd_quantile(d,rng_uniform(rng)); end do
  end subroutine ecd_random

  recursive function ecd_moment(d,order,center,lower,upper) result(v)
    type(ecd_model), intent(in) :: d
    integer, intent(in) :: order
    logical, intent(in), optional :: center
    real(dp), intent(in), optional :: lower,upper
    real(dp) :: v,c,a,b
    logical :: cen
    cen=.false.; if(present(center)) cen=center
    c=0.0_dp
    if(cen) then
      if(d%has_stats) then; c=d%stats%mean; else; c=ecd_moment(d,1,.false.); end if
    end if
    a=-huge(1.0_dp); b=huge(1.0_dp)
    if(present(lower)) a=lower
    if(present(upper)) b=upper
    v=integrate_adaptive(mfun,a,b,1e-8_dp,1e-10_dp)
  contains
    function mfun(x) result(q)
      real(dp), intent(in) :: x
      real(dp) :: q
      q=(x-c)**order*ecd_pdf(d,x)
    end function mfun
  end function ecd_moment

  function ecd_statistics(d) result(s)
    type(ecd_model), intent(in) :: d
    type(ecd_stats_type) :: s
    s%m1=ecd_moment(d,1); s%m2=ecd_moment(d,2)
    s%m3=ecd_moment(d,3); s%m4=ecd_moment(d,4)
    s%mean=s%m1; s%variance=max(0.0_dp,s%m2-s%m1*s%m1); s%stdev=sqrt(s%variance)
    if(s%variance>0.0_dp) then
      s%skewness=(s%m3-3*s%m1*s%m2+2*s%m1**3)/s%variance**1.5_dp
      s%kurtosis=(s%m4-4*s%m1*s%m3+6*s%m1*s%m1*s%m2-3*s%m1**4)/s%variance**2
    else
      s%skewness=nan_dp(); s%kurtosis=nan_dp()
    end if
  end function ecd_statistics

  function ecd_asymptotic_statistics(d,q) result(s)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: q
    type(ecd_stats_type) :: s
    real(dp) :: lo,hi,mass
    lo=ecd_quantile(d,q); hi=ecd_quantile(d,1.0_dp-q); mass=1.0_dp-2.0_dp*q
    s%m1=ecd_moment(d,1,.false.,lo,hi)/mass; s%m2=ecd_moment(d,2,.false.,lo,hi)/mass
    s%m3=ecd_moment(d,3,.false.,lo,hi)/mass; s%m4=ecd_moment(d,4,.false.,lo,hi)/mass
    s%mean=s%m1; s%variance=s%m2-s%m1*s%m1; s%stdev=sqrt(max(0.0_dp,s%variance))
    s%skewness=(s%m3-3*s%m1*s%m2+2*s%m1**3)/s%variance**1.5_dp
    s%kurtosis=(s%m4-4*s%m1*s%m3+6*s%m1*s%m1*s%m2-3*s%m1**4)/s%variance**2
  end function ecd_asymptotic_statistics

  function ecd_ellipticity(d,tol) result(r)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in), optional :: tol
    type(ellipticity_result) :: r
    real(dp) :: sd,m,dx,x,bestv,v
    integer :: i,n
    real(dp) :: tol0
    tol0=1.0e-4_dp; if(present(tol)) tol0=tol
    if(d%cusp>0 .and. d%alpha==0.0_dp .and. d%gamma==0.0_dp) return
    if(d%has_stats) then; sd=d%stats%stdev; m=d%stats%mean; else
      sd=sqrt(max(ecd_moment(d,2,.true.),tiny(1.0_dp))); m=ecd_moment(d,1)
    end if
    n=max(400,min(10000,ceiling(6.0_dp/max(tol0,1.0e-5_dp)))); dx=6.0_dp*sd/real(n,dp)
    r%left=m; bestv=-huge(1.0_dp)
    do i=0,n/2
      x=m-3.0_dp*sd+real(i,dp)*dx; v=ecd_y_slope(d,x)
      if(v>bestv) then; bestv=v; r%left=x; end if
    end do
    r%right=m; bestv=huge(1.0_dp)
    do i=n/2,n
      x=m-3.0_dp*sd+real(i,dp)*dx; v=ecd_y_slope(d,x)
      if(v<bestv) then; bestv=v; r%right=x; end if
    end do
    r%average=0.5_dp*(r%right-r%left)
  end function ecd_ellipticity

  pure elemental function ecd_cusp_std_moment(n) result(v)
    integer, intent(in) :: n
    real(dp) :: v
    if(n<0 .or. mod(n,2)/=0) then; v=0.0_dp; else; v=2.0_dp/sqrt_pi*gamma(1.5_dp*real(n+1,dp)); end if
  end function ecd_cusp_std_moment

  function ecd_cusp_std_cf(t,mu,sigma) result(v)
    real(dp), intent(in) :: t
    real(dp), intent(in), optional :: mu,sigma
    complex(dp) :: v,term
    real(dp) :: m,s
    integer :: n
    m=0.0_dp; s=1.0_dp; if(present(mu))m=mu; if(present(sigma))s=sigma
    v=(1.0_dp,0.0_dp)
    do n=2,100,2
      term=cmplx((-1.0_dp)**(n/2)*ecd_cusp_std_moment(n)*(s*t)**n/gamma(real(n+1,dp)),0.0_dp,dp)
      v=v+term
      if(abs(term)<1e-13_dp*max(1.0_dp,abs(v))) exit
    end do
    v=v*exp(cmplx(0.0_dp,m*t,dp))
  end function ecd_cusp_std_cf

  function ecd_cusp_std_mgf(t,mu,sigma) result(v)
    real(dp), intent(in) :: t
    real(dp), intent(in), optional :: mu,sigma
    real(dp) :: v,m,s,term,last
    integer :: n
    m=0.0_dp; s=1.0_dp; if(present(mu))m=mu; if(present(sigma))s=sigma
    v=1.0_dp; last=huge(1.0_dp)
    do n=2,200,2
      term=ecd_cusp_std_moment(n)*(s*t)**n/gamma(real(n+1,dp))
      if(abs(term)>last) exit
      v=v+term; last=abs(term)
      if(abs(term)<1e-13_dp*max(1.0_dp,abs(v))) exit
    end do
    v=exp(m*t)*v
  end function ecd_cusp_std_mgf

  function ecd_imgf(d,x,t,minus_one) result(v)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in), optional :: x,t
    logical, intent(in), optional :: minus_one
    real(dp) :: v,a,tt,c
    a=-huge(1.0_dp); tt=1.0_dp; c=0.0_dp
    if(present(x))a=x; if(present(t))tt=t; if(present(minus_one)) then; if(minus_one)c=1.0_dp; end if
    v=integrate_adaptive(ifun,a,huge(1.0_dp),1e-8_dp,1e-11_dp)
  contains
    function ifun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q,lg
      lg=tt*z+ecd_solve(d,z)-log(d%norm_const)
      if(lg>log(huge(1.0_dp))) then; q=huge(1.0_dp); else; q=(exp(tt*z)-c)*ecd_pdf(d,z); end if
    end function ifun
  end function ecd_imgf

  function ecd_mu_d(d) result(v)
    type(ecd_model), intent(in) :: d
    real(dp) :: v
    v=-log(ecd_imgf(d))
  end function ecd_mu_d

  function ecd_ogf(d,k,option_type) result(v)
    type(ecd_model), intent(in) :: d
    real(dp), intent(in) :: k
    character(len=*), intent(in), optional :: option_type
    real(dp) :: v
    character(len=1) :: ot
    ot='c'; if(present(option_type)) ot=option_type(1:1)
    if(ot=='c' .or. ot=='C') then
      v=integrate_adaptive(callfun,k,huge(1.0_dp),1e-8_dp,1e-11_dp)
    else
      v=integrate_adaptive(putfun,-huge(1.0_dp),k,1e-8_dp,1e-11_dp)
    end if
  contains
    function callfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=(exp(z)-exp(k))*ecd_pdf(d,z)
    end function callfun
    function putfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=(exp(k)-exp(z))*ecd_pdf(d,z)
    end function putfun
  end function ecd_ogf

  pure function ecd_y0_isomorphic(theta,radius) result(y)
    real(dp), intent(in) :: theta
    real(dp), intent(in), optional :: radius
    real(dp) :: y,r,fs,t
    r=1.0_dp; if(present(radius))r=radius
    if(abs(theta/pi-1.75_dp)<=10*epsilon(1.0_dp)) then
      y=-(r*cos(theta)/2.0_dp)**(1.0_dp/3.0_dp)
    else if(theta==0.0_dp) then
      y=r**(1.0_dp/3.0_dp)
    else if(theta>0.0_dp .and. theta<=pi) then
      fs=-2.0_dp**(2.0_dp/3.0_dp)*abs(sin(theta))**(1.0_dp/3.0_dp)*r**(1.0_dp/3.0_dp)
      t=log(tan(theta/2.0_dp))/3.0_dp; y=fs*sinh(t)
    else
      y=nan_dp()
    end if
  end function ecd_y0_isomorphic

  function ecd_max_kurtosis(j_invariant) result(ans)
    integer, intent(in), optional :: j_invariant
    real(dp) :: ans(2),parameter,kurt,best_kurt
    type(ecd_model) :: d
    type(ecd_stats_type) :: s
    integer :: jv,i,n
    jv=0; if(present(j_invariant))jv=j_invariant
    best_kurt=-huge(1.0_dp); ans=nan_dp()
    if(jv==0) then
      n=15
      do i=0,n
        parameter=2.85_dp+0.01_dp*real(i,dp)
        d=ecd_new(alpha=parameter,gamma=0.0_dp,sigma=1.0_dp)
        s=ecd_statistics(d); kurt=s%kurtosis
        if(kurt>best_kurt)then;best_kurt=kurt;ans=[parameter,kurt];end if
      end do
    else if(jv==1728) then
      n=20
      do i=0,n
        parameter=1.480_dp+0.001_dp*real(i,dp)
        d=ecd_new(alpha=0.0_dp,gamma=parameter,sigma=1.0_dp)
        s=ecd_statistics(d); kurt=s%kurtosis
        if(kurt>best_kurt)then;best_kurt=kurt;ans=[parameter,kurt];end if
      end do
    end if
  end function ecd_max_kurtosis

  function ecd_rational(x,max_den) result(frac)
    real(dp), intent(in) :: x
    integer, intent(in), optional :: max_den
    integer :: frac(2),md,p,q
    real(dp) :: err,best
    md=10000; if(present(max_den))md=max_den
    best=huge(1.0_dp); frac=[0,1]
    do q=1,md
      p=nint(x*real(q,dp)); err=abs(x-real(p,dp)/real(q,dp))
      if(err<best) then; best=err; frac=[p,q]; end if
      if(err<1e-12_dp) exit
    end do
  end function ecd_rational

end module ecd_core
