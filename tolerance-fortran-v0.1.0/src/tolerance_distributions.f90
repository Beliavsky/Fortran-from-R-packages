! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_distributions
  use tolerance_kinds, only : dp, huge_dp
  use tolerance_math, only : adaptive_simpson, beta_i, log_beta_fn, rng_beta, rng_gamma, rng_poisson, &
    hypergeom_pmf, hypergeom_cdf, clamp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private

  public :: d2exp, p2exp, q2exp, r2exp
  public :: ddpareto, pdpareto, qdpareto, rdpareto
  public :: dpoislind, ppoislind, qpoislind, rpoislind
  public :: dnhyper, pnhyper, qnhyper, rnhyper
  public :: zeta_fun, dzipfman, pzipfman, qzipfman, rzipfman
  public :: appell_f1
  public :: ddiffprop, pdiffprop, qdiffprop, rdiffprop

contains

  elemental real(dp) function d2exp(x,rate,shift,log_density) result(v)
    real(dp),intent(in)::x,rate
    real(dp),intent(in),optional::shift
    logical,intent(in),optional::log_density
    real(dp)::sh
    logical::lg
    sh=0.0_dp;if(present(shift))sh=shift
    lg=.false.;if(present(log_density))lg=log_density
    if(rate<=0.0_dp .or. x<sh)then
      v=merge(-huge_dp,0.0_dp,lg)
    else
      v=-(x-sh)/rate-log(rate)
      if(.not.lg)v=exp(v)
    end if
  end function d2exp

  elemental real(dp) function p2exp(q,rate,shift,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,rate
    real(dp),intent(in),optional::shift
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::sh
    logical::lt,lp
    sh=0.0_dp;if(present(shift))sh=shift
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    if(rate<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    if(q<sh)then;v=0.0_dp;else;v=1.0_dp-exp(-(q-sh)/rate);end if
    if(.not.lt)v=1.0_dp-v
    if(lp)v=log(v)
  end function p2exp

  elemental real(dp) function q2exp(p,rate,shift,lower_tail,log_p) result(x)
    real(dp),intent(in)::p,rate
    real(dp),intent(in),optional::shift
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp,sh
    logical::lt,lp
    sh=0.0_dp;if(present(shift))sh=shift
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1.0_dp-pp
    if(rate<=0.0_dp .or. pp<0.0_dp .or. pp>1.0_dp)then
      x=ieee_value(0.0_dp,ieee_quiet_nan)
    else if(pp>=1.0_dp)then
      x=huge_dp
    else
      x=sh-rate*log(1.0_dp-pp)
    end if
  end function q2exp

  real(dp) function r2exp(rate,shift) result(x)
    real(dp),intent(in)::rate
    real(dp),intent(in),optional::shift
    real(dp)::u,sh
    sh=0.0_dp;if(present(shift))sh=shift
    call random_number(u);x=q2exp(u,rate,sh)
  end function r2exp

  elemental real(dp) function ddpareto(x,theta,log_density) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::theta
    logical,intent(in),optional::log_density
    logical::lg
    real(dp)::a,b
    lg=.false.;if(present(log_density))lg=log_density
    if(theta<=0.0_dp .or. theta>=1.0_dp .or. x<0)then
      v=merge(-huge_dp,0.0_dp,lg);return
    end if
    a=exp(log(theta)*log(real(1+x,dp)))
    b=exp(log(theta)*log(real(2+x,dp)))
    v=max(0.0_dp,a-b)
    if(lg)then
      if(v>0.0_dp)then;v=log(v);else;v=-huge_dp;end if
    end if
  end function ddpareto

  elemental real(dp) function pdpareto(q,theta,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,theta
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    integer::k
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    if(theta<=0.0_dp .or. theta>=1.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(q<0.0_dp)then
      v=0.0_dp
    else
      k=floor(q)
      v=1.0_dp-exp(log(theta)*log(real(2+k,dp)))
    end if
    if(.not.lt)v=1.0_dp-v
    if(lp)v=log(v)
  end function pdpareto

  real(dp) function qdpareto(p,theta,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,theta
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp
    logical::lt,lp
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1.0_dp-pp
    if(theta<=0.0_dp .or. theta>=1.0_dp .or. pp<0.0_dp .or. pp>1.0_dp)then
      q=ieee_value(0.0_dp,ieee_quiet_nan)
    else if(pp>=1.0_dp)then
      q=huge_dp
    else if(pp<=0.0_dp)then
      q=0.0_dp
    else
      q=real(max(floor(exp(log(1.0_dp-pp)/log(theta))-2.0_dp),0),dp)
    end if
  end function qdpareto

  integer function rdpareto(theta) result(x)
    real(dp),intent(in)::theta
    real(dp)::u
    call random_number(u)
    x=int(qdpareto(u,theta))
  end function rdpareto

  elemental real(dp) function dpoislind(x,theta,log_density) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::theta
    logical,intent(in),optional::log_density
    logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(theta<=0.0_dp .or. x<0)then;v=merge(-huge_dp,0.0_dp,lg);return;end if
    v=2.0_dp*log(theta)+log(real(x,dp)+theta+2.0_dp)-real(x+3,dp)*log(theta+1.0_dp)
    if(.not.lg)v=exp(v)
  end function dpoislind

  real(dp) function ppoislind(q,theta,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,theta
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    integer::k,j
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    if(q<0.0_dp)then
      v=0.0_dp
    else
      k=floor(q);v=0.0_dp
      do j=0,k;v=v+dpoislind(j,theta);end do
      v=clamp(v,0.0_dp,1.0_dp)
    end if
    if(.not.lt)v=1.0_dp-v
    if(lp)v=log(v)
  end function ppoislind

  real(dp) function qpoislind(p,theta,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,theta
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp
    logical::lt,lp
    integer::lo,hi,mid
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1.0_dp-pp
    if(pp<0.0_dp .or. pp>1.0_dp .or. theta<=0.0_dp)then;q=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(pp<=0.0_dp)then;q=0.0_dp;return;end if
    if(pp>=1.0_dp)then;q=huge_dp;return;end if
    lo=0;hi=max(10,int(20.0_dp/theta+50.0_dp))
    do while(ppoislind(real(hi,dp),theta)<pp .and. hi<100000000);hi=2*hi;end do
    do while(lo<hi)
      mid=(lo+hi)/2
      if(ppoislind(real(mid,dp),theta)>=pp)then;hi=mid;else;lo=mid+1;end if
    end do
    q=real(lo,dp)
  end function qpoislind

  integer function rpoislind(theta) result(x)
    real(dp),intent(in)::theta
    real(dp)::u,lambda
    call random_number(u)
    lambda=rng_gamma(1.0_dp,1.0_dp/theta)
    if(u>theta/(theta+1.0_dp))lambda=lambda+rng_gamma(1.0_dp,1.0_dp/theta)
    x=rng_poisson(lambda)
  end function rpoislind

  real(dp) function dnhyper(x,m,n,k,log_density) result(v)
    integer,intent(in)::x,m,n,k
    logical,intent(in),optional::log_density
    logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(k>m .or. x<k .or. x>n-m+k)then;v=merge(-huge_dp,0.0_dp,lg);return;end if
    v=hypergeom_pmf(k-1,m,n-m,x-1)*real(m-k+1,dp)/real(n-x+1,dp)
    if(lg)then;if(v>0.0_dp)then;v=log(v);else;v=-huge_dp;end if;end if
  end function dnhyper

  real(dp) function pnhyper(q,m,n,k,lower_tail,log_p) result(v)
    integer,intent(in)::q,m,n,k
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lp
    integer::x
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    if(q<k)then
      v=0.0_dp
    else if(q>=n-m+k)then
      v=1.0_dp
    else
      v=0.0_dp
      do x=k,q;v=v+dnhyper(x,m,n,k);end do
      v=clamp(v,0.0_dp,1.0_dp)
    end if
    if(.not.lt)v=1.0_dp-v
    if(lp)v=log(v)
  end function pnhyper

  integer function qnhyper(p,m,n,k,lower_tail,log_p) result(q)
    real(dp),intent(in)::p
    integer,intent(in)::m,n,k
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp
    logical::lt,lp
    integer::lo,hi,mid
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1.0_dp-pp
    lo=k;hi=n-m+k
    do while(lo<hi)
      mid=(lo+hi)/2
      if(pnhyper(mid,m,n,k)>=pp)then;hi=mid;else;lo=mid+1;end if
    end do
    q=lo
  end function qnhyper

  integer function rnhyper(m,n,k) result(q)
    integer,intent(in)::m,n,k
    real(dp)::u
    call random_number(u);q=qnhyper(u,m,n,k)
  end function rnhyper

  real(dp) function zeta_fun(s,a) result(z)
    real(dp),intent(in)::s
    real(dp),intent(in),optional::a
    real(dp)::aa,x,sumv,term
    integer::k,n
    ! Euler-Maclaurin evaluation of the Hurwitz zeta; R's zeta.fun uses the
    ! Riemann zeta for a=1.  This form also supports Zipf-Mandelbrot normalizers.
    aa=1.0_dp;if(present(a))aa=a
    if(s<=1.0_dp .or. aa<=0.0_dp)then;z=huge_dp;return;end if
    n=32;sumv=0.0_dp
    do k=0,n-1;sumv=sumv/(1.0_dp) + (aa+real(k,dp))**(-s);end do
    x=aa+real(n,dp)
    z=sumv+x**(1.0_dp-s)/(s-1.0_dp)+0.5_dp*x**(-s)
    term=s*x**(-s-1.0_dp)/12.0_dp;z=z+term
    term=-s*(s+1.0_dp)*(s+2.0_dp)*x**(-s-3.0_dp)/720.0_dp;z=z+term
    term=s*(s+1.0_dp)*(s+2.0_dp)*(s+3.0_dp)*(s+4.0_dp)*x**(-s-5.0_dp)/30240.0_dp
    z=z+term
  end function zeta_fun

  real(dp) function dzipfman(x,s,b,nmax,log_density) result(v)
    integer,intent(in)::x
    real(dp),intent(in)::s
    real(dp),intent(in),optional::b
    integer,intent(in),optional::nmax
    logical,intent(in),optional::log_density
    real(dp)::bb,norm
    integer::j,nm
    logical::lg
    bb=0.0_dp;if(present(b))bb=b
    nm=-1;if(present(nmax))nm=nmax
    lg=.false.;if(present(log_density))lg=log_density
    if(x<1 .or. s<=0.0_dp .or. bb<0.0_dp .or. (nm>0 .and. x>nm))then
      v=merge(-huge_dp,0.0_dp,lg);return
    end if
    if(nm>0)then
      norm=0.0_dp;do j=1,nm;norm=norm+(real(j,dp)+bb)**(-s);end do
    else
      if(bb/=0.0_dp .or. s<=1.0_dp)then;v=merge(-huge_dp,0.0_dp,lg);return;end if
      norm=zeta_fun(s)
    end if
    v=-s*log(real(x,dp)+bb)-log(norm)
    if(.not.lg)v=exp(v)
  end function dzipfman

  real(dp) function pzipfman(q,s,b,nmax,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,s
    real(dp),intent(in),optional::b
    integer,intent(in),optional::nmax
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::bb
    integer::k,j,nm
    logical::lt,lp
    bb=0.0_dp;if(present(b))bb=b
    nm=-1;if(present(nmax))nm=nmax
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    if(q<1.0_dp)then;v=0.0_dp;else
      k=floor(q);if(nm>0)k=min(k,nm);v=0.0_dp
      do j=1,k;v=v+dzipfman(j,s,bb,nm);end do
      v=clamp(v,0.0_dp,1.0_dp)
    end if
    if(.not.lt)v=1.0_dp-v
    if(lp)v=log(v)
  end function pzipfman

  real(dp) function qzipfman(p,s,b,nmax,lower_tail,log_p) result(q)
    real(dp),intent(in)::p,s
    real(dp),intent(in),optional::b
    integer,intent(in),optional::nmax
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::pp,bb
    integer::lo,hi,mid,nm
    logical::lt,lp
    bb=0.0_dp;if(present(b))bb=b
    nm=-1;if(present(nmax))nm=nmax
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lp=.false.;if(present(log_p))lp=log_p
    pp=p;if(lp)pp=exp(pp);if(.not.lt)pp=1.0_dp-pp
    if(pp<=0.0_dp)then;q=1.0_dp;return;end if
    if(pp>=1.0_dp)then;q=merge(real(nm,dp),huge_dp,nm>0);return;end if
    lo=1;hi=merge(nm,1024,nm>0)
    if(nm<=0)then
      do while(pzipfman(real(hi,dp),s,bb,nm)<pp .and. hi<100000000);hi=2*hi;end do
    end if
    do while(lo<hi)
      mid=(lo+hi)/2
      if(pzipfman(real(mid,dp),s,bb,nm)>=pp)then;hi=mid;else;lo=mid+1;end if
    end do
    q=real(lo,dp)
  end function qzipfman

  integer function rzipfman(s,b,nmax) result(x)
    real(dp),intent(in)::s
    real(dp),intent(in),optional::b
    integer,intent(in),optional::nmax
    real(dp)::u,bb
    integer::nm
    bb=0.0_dp;if(present(b))bb=b
    nm=-1;if(present(nmax))nm=nmax
    call random_number(u);x=int(qzipfman(u,s,bb,nm))
  end function rzipfman

  real(dp) function appell_f1(a,b,bprime,c,x,y) result(v)
    real(dp),intent(in)::a,b,bprime,c,x,y
    real(dp)::coef
    if(a<=0.0_dp .or. c<=a)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    coef=exp(log_gamma(c)-log_gamma(a)-log_gamma(c-a))
    v=coef*adaptive_simpson(integrand,1.0e-12_dp,1.0_dp-1.0e-12_dp,1.0e-9_dp,24)
  contains
    function integrand(u) result(z)
      real(dp),intent(in)::u
      real(dp)::z
      z=u**(a-1.0_dp)*(1.0_dp-u)**(c-a-1.0_dp) * &
        (1.0_dp-u*x)**(-b)*(1.0_dp-u*y)**(-bprime)
    end function integrand
  end function appell_f1

  real(dp) function ddiffprop(d,k1,k2,n1,n2,a1,a2) result(v)
    real(dp),intent(in)::d
    integer,intent(in)::k1,k2,n1,n2
    real(dp),intent(in),optional::a1,a2
    real(dp)::aa1,aa2,c1,c2,b1,b2,lo,hi,lnorm
    aa1=0.5_dp;if(present(a1))aa1=a1
    aa2=0.5_dp;if(present(a2))aa2=a2
    if(d<=-1.0_dp .or. d>=1.0_dp)then;v=0.0_dp;return;end if
    c1=real(k1,dp)+aa1;b1=real(n1-k1,dp)+aa1
    c2=real(k2,dp)+aa2;b2=real(n2-k2,dp)+aa2
    lnorm=log_beta_fn(c1,b1)+log_beta_fn(c2,b2)
    lo=max(0.0_dp,-d);hi=min(1.0_dp,1.0_dp-d)
    if(hi<=lo)then;v=0.0_dp;return;end if
    v=adaptive_simpson(integrand,lo+1.0e-12_dp,hi-1.0e-12_dp,1.0e-8_dp,22)
  contains
    function integrand(y0) result(z)
      real(dp),intent(in)::y0
      real(dp)::z,x0
      x0=y0+d
      if(x0<=0.0_dp .or. x0>=1.0_dp .or. y0<=0.0_dp .or. y0>=1.0_dp)then
        z=0.0_dp
      else
        z=exp((c1-1.0_dp)*log(x0)+(b1-1.0_dp)*log(1.0_dp-x0) + &
          (c2-1.0_dp)*log(y0)+(b2-1.0_dp)*log(1.0_dp-y0)-lnorm)
      end if
    end function integrand
  end function ddiffprop

  real(dp) function pdiffprop(q,k1,k2,n1,n2,a1,a2) result(v)
    real(dp),intent(in)::q
    integer,intent(in)::k1,k2,n1,n2
    real(dp),intent(in),optional::a1,a2
    real(dp)::aa1,aa2
    aa1=0.5_dp;if(present(a1))aa1=a1
    aa2=0.5_dp;if(present(a2))aa2=a2
    if(q<=-1.0_dp)then;v=0.0_dp;return;end if
    if(q>=1.0_dp)then;v=1.0_dp;return;end if
    v=adaptive_simpson(integrand,-1.0_dp+1.0e-9_dp,q,2.0e-7_dp,20)
    v=clamp(v,0.0_dp,1.0_dp)
  contains
    function integrand(d) result(z)
      real(dp),intent(in)::d
      real(dp)::z
      z=ddiffprop(d,k1,k2,n1,n2,aa1,aa2)
    end function integrand
  end function pdiffprop

  real(dp) function qdiffprop(p,k1,k2,n1,n2,a1,a2) result(q)
    real(dp),intent(in)::p
    integer,intent(in)::k1,k2,n1,n2
    real(dp),intent(in),optional::a1,a2
    real(dp)::aa1,aa2,lo,hi,mid
    integer::i
    aa1=0.5_dp;if(present(a1))aa1=a1
    aa2=0.5_dp;if(present(a2))aa2=a2
    if(p<=0.0_dp)then;q=-1.0_dp;return;end if
    if(p>=1.0_dp)then;q=1.0_dp;return;end if
    lo=-1.0_dp;hi=1.0_dp
    do i=1,70
      mid=0.5_dp*(lo+hi)
      if(pdiffprop(mid,k1,k2,n1,n2,aa1,aa2)<p)then;lo=mid;else;hi=mid;end if
    end do
    q=0.5_dp*(lo+hi)
  end function qdiffprop

  real(dp) function rdiffprop(k1,k2,n1,n2,a1,a2) result(x)
    integer,intent(in)::k1,k2,n1,n2
    real(dp),intent(in),optional::a1,a2
    real(dp)::aa1,aa2
    aa1=0.5_dp;if(present(a1))aa1=a1
    aa2=0.5_dp;if(present(a2))aa2=a2
    x=rng_beta(real(k1,dp)+aa1,real(n1-k1,dp)+aa1) - &
      rng_beta(real(k2,dp)+aa2,real(n2-k2,dp)+aa2)
  end function rdiffprop

end module tolerance_distributions
