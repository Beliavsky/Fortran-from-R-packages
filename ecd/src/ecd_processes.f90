! SPDX-License-Identifier: Artistic-2.0
module ecd_processes
  use ecd_kinds, only : dp, pi, sqrt_pi, ecd_ok, ecd_invalid
  use ecd_math, only : nan_dp, integrate_adaptive, brent_root, bessel_k, gamma_p, gamma_quantile
  use ecd_rng, only : rng_state, rng_uniform, rng_normal, rng_exponential, rng_gamma
  implicit none
  private

  type, public :: cumulants4
    real(dp) :: mean=0.0_dp, variance=0.0_dp, k3=0.0_dp, k4=0.0_dp
    real(dp) :: skewness=0.0_dp, kurtosis=3.0_dp
  end type cumulants4

  type, public :: sld_model
    real(dp) :: t=1.0_dp, nu0=0.0_dp, theta=1.0_dp, convo=1.0_dp
    real(dp) :: beta_a=0.0_dp, mu=0.0_dp, lambda=4.0_dp
  end type sld_model

  public :: laplace_pdf, laplace_random
  public :: stdlap_pdf, stdlap_cdf, stdlap_quantile, stdlap_random, stdlap_cf, stdlap_cumulants
  public :: stdlap_pdf_poly
  public :: stable_count_pdf, stable_count_cdf, stable_count_quantile, stable_count_random
  public :: stable_count_cf, stable_count_cumulants
  public :: sld_new, sld_pdf, sld_cdf, sld_quantile, sld_random, sld_cf, sld_cumulants
  public :: qsl_variance_analytic, qsl_skewness_analytic, qsl_kurtosis_analytic
  public :: qsl_std_pdf0_analytic, qsl_pdf_integrand_analytic
  public :: levy_dlambda, levy_dskewed, k2moments, moments2k
  public :: stable_pdf_positive

contains

  pure elemental function laplace_pdf(x,b) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b
    real(dp) :: v,sc
    sc=1.0_dp; if(present(b))sc=b
    if(sc<=0.0_dp) then; v=nan_dp(); else; v=0.5_dp/sc*exp(-abs(x)/sc); end if
  end function laplace_pdf

  subroutine laplace_random(rng,x,b)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: b
    real(dp) :: sc
    integer :: i
    sc=1.0_dp; if(present(b))sc=b
    do i=1,size(x)
      x(i)=rng_exponential(rng,sc)*merge(1.0_dp,-1.0_dp,rng_uniform(rng)>=0.5_dp)
    end do
  end subroutine laplace_random

  function stdlap_pdf(x,t,convo,beta,mu) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: t,convo,beta,mu
    real(dp) :: v,tt,m,b,loc,sm,b0,w,c,k,c0,ck,eb,z
    tt=1.0_dp; m=1.0_dp; b=0.0_dp; loc=0.0_dp
    if(present(t))tt=t; if(present(convo))m=convo; if(present(beta))b=beta; if(present(mu))loc=mu
    if(tt<=0.0_dp .or. m<=0.0_dp) then; v=nan_dp(); return; end if
    z=x-loc; sm=sqrt(tt/(2.0_dp+b*b)/m); b0=sqrt(1.0_dp+b*b/4.0_dp); w=abs(b0*z/sm)
    c=1.0_dp/gamma(m)/sqrt_pi/sm
    if(w==0.0_dp) then
      if(m>0.5_dp) then; c0=c*gamma(m-0.5_dp)/2.0_dp; ck=c0
      else; v=huge(1.0_dp); return; end if
    else
      k=(w/2.0_dp)**(m-0.5_dp)*bessel_k(m-0.5_dp,w); ck=c*k
    end if
    eb=(b0**(-2.0_dp))**(m-0.5_dp)*exp(b*z/(2.0_dp*sm))
    v=ck*eb
  end function stdlap_pdf

  function stdlap_cdf(x,t,convo,beta,mu) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: t,convo,beta,mu
    real(dp) :: v,loc
    loc=0.0_dp; if(present(mu))loc=mu
    if(x<=loc) then
      v=integrate_adaptive(pdf_fun,-huge(1.0_dp),x,1e-9_dp,1e-11_dp)
    else
      v=1.0_dp-integrate_adaptive(pdf_fun,x,huge(1.0_dp),1e-9_dp,1e-11_dp)
    end if
    v=max(0.0_dp,min(1.0_dp,v))
  contains
    function pdf_fun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=stdlap_pdf(z,t,convo,beta,mu)
    end function pdf_fun
  end function stdlap_cdf

  function stdlap_quantile(p,t,convo,beta,mu,status) result(x)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: t,convo,beta,mu
    integer, intent(out), optional :: status
    real(dp) :: x,tt,loc,lo,hi
    integer :: st
    tt=1.0_dp; loc=0.0_dp; if(present(t))tt=t; if(present(mu))loc=mu
    if(present(status))status=ecd_ok
    if(p<0.0_dp .or. p>1.0_dp) then; x=nan_dp(); if(present(status))status=ecd_invalid; return; end if
    if(p==0.0_dp) then; x=-huge(1.0_dp); return; end if
    if(p==1.0_dp) then; x=huge(1.0_dp); return; end if
    lo=loc-20.0_dp*sqrt(tt); hi=loc+20.0_dp*sqrt(tt)
    x=brent_root(qfun,lo,hi,1e-10_dp,200,st)
    if(present(status))status=st
  contains
    function qfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=stdlap_cdf(z,t,convo,beta,mu)-p
    end function qfun
  end function stdlap_quantile

  subroutine stdlap_random(rng,x,t,convo,beta,mu)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: t,convo,beta,mu
    real(dp) :: tt,m,b,loc,sm,g
    integer :: i
    tt=1.0_dp; m=1.0_dp; b=0.0_dp; loc=0.0_dp
    if(present(t))tt=t; if(present(convo))m=convo; if(present(beta))b=beta; if(present(mu))loc=mu
    sm=sqrt(tt/(2.0_dp+b*b)/m)
    do i=1,size(x)
      g=rng_gamma(rng,m,1.0_dp)
      x(i)=loc+b*sm*g+sqrt(2.0_dp)*sm*sqrt(g)*rng_normal(rng)
    end do
  end subroutine stdlap_random

  pure function stdlap_cf(s,t,convo,beta,mu) result(v)
    real(dp), intent(in) :: s
    real(dp), intent(in), optional :: t,convo,beta,mu
    complex(dp) :: v,base
    real(dp) :: tt,m,b,loc,sm
    tt=1.0_dp; m=1.0_dp; b=0.0_dp; loc=0.0_dp
    if(present(t))tt=t; if(present(convo))m=convo; if(present(beta))b=beta; if(present(mu))loc=mu
    sm=sqrt(tt/(2.0_dp+b*b)/m)
    base=cmplx(1.0_dp+s*s*sm*sm,-b*s*sm,dp)
    v=exp(cmplx(0.0_dp,loc*s,dp))*base**(-m)
  end function stdlap_cf

  pure function stdlap_cumulants(t,convo,beta,mu) result(k)
    real(dp), intent(in), optional :: t,convo,beta,mu
    type(cumulants4) :: k
    real(dp) :: tt,m,b,loc,bb
    tt=1.0_dp; m=1.0_dp; b=0.0_dp; loc=0.0_dp
    if(present(t))tt=t; if(present(convo))m=convo; if(present(beta))b=beta; if(present(mu))loc=mu
    bb=2.0_dp+b*b
    k%mean=loc+b*sqrt(m*tt/bb)
    k%variance=tt
    k%k3=2.0_dp*tt**1.5_dp/sqrt(m)*b*(3.0_dp+b*b)/bb**1.5_dp
    k%k4=6.0_dp*tt*tt/m*(2.0_dp+4.0_dp*b*b+b**4)/bb**2
    k%skewness=k%k3/k%variance**1.5_dp; k%kurtosis=k%k4/k%variance**2+3.0_dp
  end function stdlap_cumulants

  function stdlap_pdf_poly(x,t,convo,beta,mu) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: t,convo,beta,mu
    real(dp) :: v,tt,m,b,loc,z,sm,b0,bp,bn,w,bs,ey,p
    integer :: mi
    tt=1.0_dp; m=1.0_dp; b=0.0_dp; loc=0.0_dp
    if(present(t))tt=t; if(present(convo))m=convo; if(present(beta))b=beta; if(present(mu))loc=mu
    mi=nint(m)
    if(abs(m-real(mi,dp))>1e-12_dp .or. mi<1 .or. mi>4) then; v=nan_dp(); return; end if
    z=x-loc; sm=sqrt(tt/(2+b*b)/m); b0=sqrt(1+b*b/4); bp=b0+b/2; bn=b0-b/2
    w=abs(b0*z/sm); bs=merge(bp,bn,z>0); ey=exp(-abs(z/sm/bs))
    select case(mi)
    case(1); v=ey/sm/(2*b0)
    case(2); v=ey/sm*(w+1)/(4*b0**3)
    case(3); p=w*w+3*w+3; v=ey/sm*p/(16*b0**5)
    case(4); p=w**3+6*w*w+15*w+15; v=ey/sm*p/(96*b0**7)
    end select
  end function stdlap_pdf_poly

  function stable_pdf_positive(x,alpha) result(v)
    real(dp), intent(in) :: x,alpha
    real(dp) :: v,c,s
    if(x<=0.0_dp .or. alpha<=0.0_dp .or. alpha>=1.0_dp) then; v=0.0_dp; return; end if
    c=cos(pi*alpha/2.0_dp); s=sin(pi*alpha/2.0_dp)
    v=integrate_adaptive(fourier,0.0_dp,upper_t(),2e-7_dp,1e-10_dp)/pi
    if(v<0.0_dp .and. abs(v)<1e-10_dp)v=0.0_dp
  contains
    function fourier(t) result(q)
      real(dp), intent(in) :: t
      real(dp) :: q,ta
      ta=t**alpha
      q=exp(-c*ta)*cos(x*t-s*ta)
    end function fourier
    function upper_t() result(u)
      real(dp) :: u
      u=max(40.0_dp,(35.0_dp/max(c,1e-6_dp))**(1.0_dp/alpha))
      u=min(u,2.0e5_dp)
    end function upper_t
  end function stable_pdf_positive

  function stable_count_pdf(x,alpha,nu0,theta,lambda) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha,nu0,theta,lambda
    real(dp) :: v,a,n0,th,x0,c
    a=0.5_dp; n0=0.0_dp; th=1.0_dp
    if(present(alpha))a=alpha
    if(present(lambda))a=2.0_dp/lambda
    if(present(nu0))n0=nu0
    if(present(theta))th=theta
    if(a<=0.0_dp .or. a>=1.0_dp .or. th<=0.0_dp) then; v=nan_dp(); return; end if
    x0=x-n0
    if(x0<0.0_dp) then; v=0.0_dp; return; end if
    if(abs(a-0.5_dp)<1e-14_dp) then
      if(x0==0.0_dp) then; v=0.0_dp
      else; v=x0**0.5_dp*exp(-x0/(4*th))/(4*sqrt_pi*th**1.5_dp); end if
    else
      if(x0==0.0_dp) then; v=0.0_dp; return; end if
      c=a/gamma(1.0_dp/a)
      v=c/x0*stable_pdf_positive(th/x0,a)
    end if
  end function stable_count_pdf

  function stable_count_cdf(x,alpha,nu0,theta,lambda) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha,nu0,theta,lambda
    real(dp) :: v,a,n0,th
    a=0.5_dp; n0=0.0_dp; th=1.0_dp
    if(present(alpha))a=alpha; if(present(lambda))a=2.0_dp/lambda
    if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(x<=n0) then; v=0.0_dp
    else if(abs(a-0.5_dp)<1e-14_dp) then; v=gamma_p(1.5_dp,(x-n0)/(4*th))
    else; v=integrate_adaptive(pdf_fun,n0,x,2e-6_dp,1e-9_dp); end if
    v=max(0.0_dp,min(1.0_dp,v))
  contains
    function pdf_fun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=stable_count_pdf(z,a,n0,th)
    end function pdf_fun
  end function stable_count_cdf

  function stable_count_quantile(p,alpha,nu0,theta,lambda,status) result(x)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha,nu0,theta,lambda
    integer, intent(out), optional :: status
    real(dp) :: x,a,n0,th,lo,hi
    integer :: st
    a=0.5_dp; n0=0.0_dp; th=1.0_dp
    if(present(alpha))a=alpha; if(present(lambda))a=2.0_dp/lambda
    if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(status))status=ecd_ok
    if(p<0.0_dp .or. p>1.0_dp) then; x=nan_dp(); if(present(status))status=ecd_invalid; return; end if
    if(p==0.0_dp) then; x=n0; return; end if
    if(p==1.0_dp) then; x=huge(1.0_dp); return; end if
    if(abs(a-0.5_dp)<1e-14_dp) then
      x=n0+gamma_quantile(p,1.5_dp,4*th,st); if(present(status))status=st; return
    end if
    lo=n0+1e-10_dp*th; hi=n0+100*th
    do while(stable_count_cdf(hi,a,n0,th)<p .and. hi<n0+1e8_dp*th); hi=n0+2*(hi-n0); end do
    x=brent_root(qfun,lo,hi,1e-8_dp,200,st); if(present(status))status=st
  contains
    function qfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=stable_count_cdf(z,a,n0,th)-p
    end function qfun
  end function stable_count_quantile

  subroutine stable_count_random(rng,x,alpha,nu0,theta,lambda)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: alpha,nu0,theta,lambda
    real(dp) :: a,n0,th
    integer :: i
    a=0.5_dp; n0=0.0_dp; th=1.0_dp
    if(present(alpha))a=alpha; if(present(lambda))a=2.0_dp/lambda
    if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(abs(a-0.5_dp)<1e-14_dp) then
      do i=1,size(x); x(i)=n0+rng_gamma(rng,1.5_dp,4*th); end do
    else
      do i=1,size(x); x(i)=stable_count_quantile(rng_uniform(rng),a,n0,th); end do
    end if
  end subroutine stable_count_random

  pure function stable_count_cf(s,alpha,nu0,theta,lambda) result(v)
    real(dp), intent(in) :: s
    real(dp), intent(in), optional :: alpha,nu0,theta,lambda
    complex(dp) :: v
    real(dp) :: a,n0,th
    a=0.5_dp; n0=0.0_dp; th=1.0_dp
    if(present(alpha))a=alpha; if(present(lambda))a=2.0_dp/lambda
    if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(abs(a-0.5_dp)<1e-14_dp) then
      v=exp(cmplx(0.0_dp,s*n0,dp))*cmplx(1.0_dp,-4*s*th,dp)**(-1.5_dp)
    else
      v=cmplx(nan_dp(),nan_dp(),dp)
    end if
  end function stable_count_cf

  function stable_count_cumulants(alpha,nu0,theta,lambda,status) result(k)
    real(dp), intent(in), optional :: alpha,nu0,theta,lambda
    integer, intent(out), optional :: status
    type(cumulants4) :: k
    real(dp) :: a,n0,th
    a=0.5_dp; n0=0.0_dp; th=1.0_dp
    if(present(alpha))a=alpha; if(present(lambda))a=2.0_dp/lambda
    if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(status))status=ecd_ok
    if(abs(a-0.5_dp)>1e-14_dp) then
      k%mean=nan_dp(); k%variance=nan_dp(); k%k3=nan_dp(); k%k4=nan_dp()
      k%skewness=nan_dp(); k%kurtosis=nan_dp(); if(present(status))status=ecd_invalid; return
    end if
    k%mean=n0+6*th; k%variance=24*th*th; k%k3=192*th**3; k%k4=2304*th**4
    k%skewness=k%k3/k%variance**1.5_dp; k%kurtosis=k%k4/k%variance**2+3
  end function stable_count_cumulants

  pure function sld_new(t,nu0,theta,convo,beta_a,mu,lambda) result(d)
    real(dp), intent(in), optional :: t,nu0,theta,convo,beta_a,mu,lambda
    type(sld_model) :: d
    if(present(t))d%t=t; if(present(nu0))d%nu0=nu0; if(present(theta))d%theta=theta
    if(present(convo))d%convo=convo; if(present(beta_a))d%beta_a=beta_a
    if(present(mu))d%mu=mu; if(present(lambda))d%lambda=lambda
  end function sld_new

  function sld_pdf(d,x) result(v)
    type(sld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v
    v=integrate_adaptive(nfun,0.0_dp,huge(1.0_dp),2e-7_dp,1e-10_dp)
  contains
    function nfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q,den,nu
      nu=d%theta*z; den=nu+d%nu0
      if(den<=0.0_dp) then; q=0.0_dp; return; end if
      q=d%theta*stable_count_pdf(nu,lambda=d%lambda,nu0=0.0_dp,theta=d%theta)* &
        stdlap_pdf((x-d%mu)/den,d%t,d%convo,d%beta_a*sqrt(d%t))/den
      if(.not.(q>=0.0_dp))q=0.0_dp
    end function nfun
  end function sld_pdf

  function sld_cdf(d,x) result(v)
    type(sld_model), intent(in) :: d
    real(dp), intent(in) :: x
    real(dp) :: v
    v=integrate_adaptive(nfun,0.0_dp,huge(1.0_dp),2e-7_dp,1e-9_dp)
    v=max(0.0_dp,min(1.0_dp,v))
  contains
    function nfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q,den,nu
      nu=d%theta*z; den=nu+d%nu0
      if(den<=0.0_dp) then; q=0.0_dp; return; end if
      q=d%theta*stable_count_pdf(nu,lambda=d%lambda,nu0=0.0_dp,theta=d%theta)* &
        stdlap_cdf((x-d%mu)/den,d%t,d%convo,d%beta_a*sqrt(d%t))
    end function nfun
  end function sld_cdf

  function sld_quantile(d,p,status) result(x)
    type(sld_model), intent(in) :: d
    real(dp), intent(in) :: p
    integer, intent(out), optional :: status
    real(dp) :: x,sd,lo,hi
    integer :: st
    type(cumulants4) :: ck
    ck=sld_cumulants(d)
    sd=sqrt(max(ck%variance,tiny(1.0_dp)))
    lo=d%mu-50*sd; hi=d%mu+50*sd
    x=brent_root(qfun,lo,hi,1e-7_dp,100,st)
    if(present(status))status=st
  contains
    function qfun(z) result(q)
      real(dp), intent(in) :: z
      real(dp) :: q
      q=sld_cdf(d,z)-p
    end function qfun
  end function sld_quantile

  subroutine sld_random(d,rng,x)
    type(sld_model), intent(in) :: d
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    real(dp), allocatable :: l(:),n(:)
    allocate(l(size(x)),n(size(x)))
    call stdlap_random(rng,l,d%t,d%convo,d%beta_a*sqrt(d%t))
    call stable_count_random(rng,n,lambda=d%lambda,nu0=d%nu0,theta=d%theta)
    x=l*n+d%mu
  end subroutine sld_random

  function sld_cf(d,s) result(v)
    type(sld_model), intent(in) :: d
    real(dp), intent(in) :: s
    complex(dp) :: v
    if(s==0.0_dp) then; v=cmplx(1.0_dp,0.0_dp,dp); return; end if
    v=exp(cmplx(0.0_dp,d%mu*s,dp))*integrate_adaptive_complex(cfint,0.0_dp,huge(1.0_dp))
  contains
    function cfint(z) result(q)
      real(dp), intent(in) :: z
      complex(dp) :: q
      real(dp) :: nu
      nu=d%theta*z
      q=d%theta*stdlap_cf(s*(nu+d%nu0),d%t,d%convo,d%beta_a*sqrt(d%t))* &
        stable_count_pdf(nu,lambda=d%lambda,nu0=0.0_dp,theta=d%theta)
    end function cfint
  end function sld_cf

  function integrate_adaptive_complex(f,a,b) result(v)
    interface
      function f(x) result(y)
        import dp
        real(dp), intent(in) :: x
        complex(dp) :: y
      end function f
    end interface
    real(dp), intent(in) :: a,b
    complex(dp) :: v
    v=cmplx(integrate_adaptive(fr,a,b,2e-6_dp,1e-10_dp), &
      integrate_adaptive(fi,a,b,2e-6_dp,1e-10_dp),dp)
  contains
    function fr(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y=real(f(x),dp)
    end function fr
    function fi(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y=aimag(f(x))
    end function fi
  end function integrate_adaptive_complex

  pure function k2moments(k) result(m)
    real(dp), intent(in) :: k(4)
    real(dp) :: m(4)
    m(1)=k(1); m(2)=k(2)+k(1)**2
    m(3)=k(3)+3*k(2)*k(1)+k(1)**3
    m(4)=k(4)+4*k(3)*k(1)+3*k(2)**2+6*k(2)*k(1)**2+k(1)**4
  end function k2moments

  pure function moments2k(m) result(k)
    real(dp), intent(in) :: m(4)
    real(dp) :: k(4)
    k(1)=m(1); k(2)=m(2)-m(1)**2
    k(3)=m(3)-3*m(2)*m(1)+2*m(1)**3
    k(4)=m(4)-4*m(3)*m(1)-3*m(2)**2+12*m(2)*m(1)**2-6*m(1)**4
  end function moments2k

  function sld_cumulants(d,status) result(k)
    type(sld_model), intent(in) :: d
    integer, intent(out), optional :: status
    type(cumulants4) :: k,kl,kn
    real(dp) :: ml(4),mn(4),mc(4),kc(4)
    integer :: st
    kl=stdlap_cumulants(d%t,d%convo,d%beta_a*sqrt(d%t))
    kn=stable_count_cumulants(lambda=d%lambda,nu0=d%nu0,theta=d%theta,status=st)
    if(present(status))status=st
    if(st/=ecd_ok) then
      k%mean=nan_dp(); k%variance=nan_dp(); k%k3=nan_dp(); k%k4=nan_dp();
      k%skewness=nan_dp(); k%kurtosis=nan_dp(); return
    end if
    ml=k2moments([kl%mean,kl%variance,kl%k3,kl%k4])
    mn=k2moments([kn%mean,kn%variance,kn%k3,kn%k4])
    mc=ml*mn; kc=moments2k(mc)
    k%mean=kc(1)+d%mu; k%variance=kc(2); k%k3=kc(3); k%k4=kc(4)
    k%skewness=kc(3)/kc(2)**1.5_dp; k%kurtosis=kc(4)/kc(2)**2+3
  end function sld_cumulants

  pure function qsl_variance_analytic(t,nu0,theta,convo,beta_a) result(v)
    real(dp), intent(in), optional :: t,nu0,theta,convo,beta_a
    real(dp) :: v,tt,n0,th,m,b,nu1
    tt=1.; n0=0.; th=1.; m=1.; b=0.
    if(present(t))tt=t; if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(convo))m=convo; if(present(beta_a))b=beta_a*sqrt(tt)
    nu1=n0+6*th
    v=((nu1**2+24*th**2)+m*b*b/(2+b*b)*(24*th**2))*tt
  end function qsl_variance_analytic

  pure function qsl_skewness_analytic(t,nu0,theta,convo,beta_a) result(v)
    real(dp), intent(in), optional :: t,nu0,theta,convo,beta_a
    real(dp) :: v,tt,n0,th,m,b,nu1,p,q,a,bb
    tt=1.; n0=0.; th=1.; m=1.; b=0.
    if(present(t))tt=t; if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(convo))m=convo; if(present(beta_a))b=beta_a*sqrt(tt)
    nu1=n0+6*th
    p=144*b*b+432+144*m*(2+b*b)
    q=192*b*b*m*m+576*(2+b*b)*m+384*b*b+1152
    a=2*(3+b*b)*nu1**3+p*th**2*nu1+q*th**3
    bb=(2+b*b)*nu1**2+24*th**2*(2+b*b*(1+m))
    v=b/sqrt(m)*a/bb**1.5_dp
  end function qsl_skewness_analytic

  pure function qsl_kurtosis_analytic(t,nu0,theta,convo,beta_a) result(v)
    real(dp), intent(in), optional :: t,nu0,theta,convo,beta_a
    real(dp) :: v,tt,n0,th,m,b,nu1,a,bb
    tt=1.; n0=0.; th=1.; m=1.; b=0.
    if(present(t))tt=t; if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(convo))m=convo; if(present(beta_a))b=beta_a
    if(b/=0.0_dp) then; v=nan_dp(); return; end if
    nu1=n0+6*th
    a=nu1**4+(144+96*m)*th**2*nu1**2+768*(1+m)*th**3*nu1+(4032+3456*m)*th**4
    bb=(nu1**2+24*th**2)**2
    v=3+(3/m)*a/bb
  end function qsl_kurtosis_analytic

  function qsl_std_pdf0_analytic(t,nu0,theta,convo,beta_a) result(v)
    real(dp), intent(in), optional :: t,nu0,theta,convo,beta_a
    real(dp) :: v,tt,n0,th,m,b,c,n,pdf0,varann
    tt=1.; n0=0.; th=1.; m=1.; b=0.
    if(present(t))tt=t; if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(convo))m=convo; if(present(beta_a))b=beta_a
    if(b/=0.0_dp) then; v=nan_dp(); return; end if
    c=gamma(m-0.5_dp)*sqrt(2*m)/gamma(m)/8/pi
    n=n0/th; pdf0=2*sqrt_pi-pi*sqrt(n)*exp(n/4)*erfc(sqrt(n/4))
    varann=(n+6)**2+24
    v=c*pdf0*sqrt(varann)
  end function qsl_std_pdf0_analytic

  function qsl_pdf_integrand_analytic(x,nu,t,nu0,theta,convo,beta_a,mu) result(v)
    real(dp), intent(in) :: x,nu
    real(dp), intent(in), optional :: t,nu0,theta,convo,beta_a,mu
    real(dp) :: v,tt,n0,th,m,b,loc,b0,sm,w,c,a,k,bb
    tt=1.; n0=0.; th=1.; m=1.; b=0.; loc=0.
    if(present(t))tt=t; if(present(nu0))n0=nu0; if(present(theta))th=theta
    if(present(convo))m=convo; if(present(beta_a))b=beta_a*sqrt(tt); if(present(mu))loc=mu
    b0=sqrt(1+b*b/4); sm=sqrt(tt/m/(2+b*b)); w=b0*(x-loc)/sm/(nu+n0)
    c=4*pi*gamma(m)*sm*th**1.5_dp
    a=sqrt(max(nu,0.0_dp))/(nu+n0)
    if(w==0.0_dp) then; k=0.0_dp; else; k=abs(w/(2*b0*b0))**(m-0.5_dp)*bessel_k(m-0.5_dp,abs(w)); end if
    bb=exp(w*b/(2*b0)-nu/(4*th)); v=a*k*bb/c
  end function qsl_pdf_integrand_analytic

  pure elemental function levy_dlambda(x,lambda) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda
    real(dp) :: v,l
    l=4.0_dp; if(present(lambda))l=lambda
    v=exp(-(x*x)**(1.0_dp/l))/(l*gamma(l/2.0_dp))
  end function levy_dlambda

  pure elemental function levy_dskewed(x,lambda) result(v)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda
    real(dp) :: v,l
    l=4.0_dp; if(present(lambda))l=lambda
    if(l==4.0_dp .and. x>0.0_dp) then
      v=exp(-1.0_dp/(4*x))/(2*sqrt_pi*x**1.5_dp)
    else
      v=0.0_dp
    end if
  end function levy_dskewed

end module ecd_processes
