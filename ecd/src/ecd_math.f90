! SPDX-License-Identifier: Artistic-2.0
module ecd_math
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use ecd_kinds, only : dp, pi, sqrt_pi, ecd_ok, ecd_invalid, ecd_no_convergence
  implicit none
  private

  public :: nan_dp, cbrt_real, normal_cdf, normal_quantile
  public :: gamma_p, gamma_q, gamma_quantile, upper_incomplete_gamma
  public :: integrate_adaptive, brent_root, bessel_k, dawson_f, erfi_f, erfcx_f
  public :: log_gamma_sign, rational_approx, hypergeom_2f0, erfq, erfq_sum

  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
  end interface

contains

  pure function nan_dp() result(x)
    real(dp) :: x
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure function cbrt_real(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if (x == 0.0_dp) then
      y = 0.0_dp
    else
      y = sign(abs(x)**(1.0_dp/3.0_dp), x)
    end if
  end function cbrt_real

  pure function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, -2.759285104469687e2_dp, &
       1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, -1.556989798598866e2_dp, &
       6.680131188771972e1_dp, -1.328068155288572e1_dp]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, -2.400758277161838_dp, &
      -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, 2.445134137142996_dp, &
       3.754408661907416_dp]
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  function gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p, ap, del, sumv, b, c, d, h, an
    integer :: n
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = nan_dp(); return
    end if
    if (x == 0.0_dp) then
      p = 0.0_dp; return
    end if
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n=1,10000
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n=1,10000
        an = -real(n,dp)*(real(n,dp)-a)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      p = 1.0_dp - exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p = max(0.0_dp,min(1.0_dp,p))
  end function gamma_p

  function gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    real(dp) :: q
    q = 1.0_dp-gamma_p(a,x)
  end function gamma_q

  function upper_incomplete_gamma(a,x) result(g)
    real(dp), intent(in) :: a,x
    real(dp) :: g
    g = gamma_q(a,x)*gamma(a)
  end function upper_incomplete_gamma

  function gamma_quantile(p, shape, scale, status) result(x)
    real(dp), intent(in) :: p, shape, scale
    integer, intent(out), optional :: status
    real(dp) :: x, lo, hi, f
    integer :: i
    if (present(status)) status = ecd_ok
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp(); if (present(status)) status=ecd_invalid; return
    end if
    if (p == 0.0_dp) then; x=0.0_dp; return; end if
    if (p == 1.0_dp) then; x=huge(1.0_dp); return; end if
    lo = 0.0_dp
    hi = max(1.0_dp, shape*scale + 10.0_dp*sqrt(shape)*scale)
    do while (gamma_p(shape,hi/scale) < p)
      hi = hi*2.0_dp
      if (hi > 1.0e300_dp) exit
    end do
    do i=1,200
      x = 0.5_dp*(lo+hi)
      f = gamma_p(shape,x/scale)
      if (f < p) then; lo=x; else; hi=x; end if
      if (abs(hi-lo) <= 1.0e-12_dp*(1.0_dp+abs(x))) exit
    end do
    x = 0.5_dp*(lo+hi)
  end function gamma_quantile

  recursive function integrate_adaptive(f, a, b, rel_tol, abs_tol, status) result(value)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: rel_tol, abs_tol
    integer, intent(out), optional :: status
    real(dp) :: value, rt, at, aa, bb
    logical :: a_finite, b_finite
    rt = 1.0e-9_dp; at = 1.0e-11_dp
    if (present(rel_tol)) rt = rel_tol
    if (present(abs_tol)) at = abs_tol
    if (present(status)) status = ecd_ok
    if (a == b) then; value=0.0_dp; return; end if
    ! Many callers use +/-huge() as portable infinity sentinels.  Treat
    ! values near the representable limit as infinite before squaring or
    ! otherwise transforming them in an integrand.
    a_finite = ieee_is_finite(a) .and. abs(a) < 0.5_dp*huge(1.0_dp)
    b_finite = ieee_is_finite(b) .and. abs(b) < 0.5_dp*huge(1.0_dp)
    if (a_finite .and. b_finite) then
      value = simpson_driver(f,a,b,rt,at)
    else if (.not.a_finite .and. .not.b_finite) then
      aa=1.0e-10_dp; bb=1.0_dp-1.0e-10_dp
      value = simpson_driver(both_inf,aa,bb,rt,at)
    else if (a_finite) then
      aa=0.0_dp; bb=1.0_dp-1.0e-10_dp
      value = simpson_driver(right_inf,aa,bb,rt,at)
    else
      aa=0.0_dp; bb=1.0_dp-1.0e-10_dp
      value = simpson_driver(left_inf,aa,bb,rt,at)
    end if
  contains
    function both_inf(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: y,x,c
      c=cos(pi*(t-0.5_dp)); x=tan(pi*(t-0.5_dp))
      y=f(x)*pi/(c*c)
    end function both_inf
    function right_inf(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: y,x
      x=a+t/(1.0_dp-t)
      y=f(x)/(1.0_dp-t)**2
    end function right_inf
    function left_inf(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: y,x
      x=b-t/(1.0_dp-t)
      y=f(x)/(1.0_dp-t)**2
    end function left_inf
  end function integrate_adaptive

  function simpson_driver(f,a,b,rt,at) result(v)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a,b,rt,at
    real(dp) :: v,fa,fb,fc,s
    fa=f(a); fb=f(b); fc=f(0.5_dp*(a+b))
    s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    v=adapt(f,a,b,fa,fb,fc,s,max(at,rt*abs(s)),24)
  end function simpson_driver

  recursive function adapt(f,a,b,fa,fb,fc,s,tol,depth) result(v)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: a,b,fa,fb,fc,s,tol
    integer, intent(in) :: depth
    real(dp) :: v,c,d,e,fd,fe,sl,sr,s2
    c=0.5_dp*(a+b); d=0.5_dp*(a+c); e=0.5_dp*(c+b)
    fd=f(d); fe=f(e)
    sl=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
    sr=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
    s2=sl+sr
    if (depth<=0 .or. abs(s2-s)<=15.0_dp*tol) then
      v=s2+(s2-s)/15.0_dp
    else
      v=adapt(f,a,c,fa,fc,fd,sl,tol/2.0_dp,depth-1)+ &
        adapt(f,c,b,fc,fb,fe,sr,tol/2.0_dp,depth-1)
    end if
  end function adapt

  recursive function brent_root(f, ax, bx, tol, max_iter, status) result(root)
    procedure(scalar_function) :: f
    real(dp), intent(in) :: ax,bx
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    integer, intent(out), optional :: status
    real(dp) :: root,lo,hi,flo,fhi,fc,c,t,x_tol
    integer :: it,nmax
    t=1.0e-10_dp; nmax=200
    if(present(tol))t=tol
    if(present(max_iter))nmax=max_iter
    if(present(status))status=ecd_ok
    lo=min(ax,bx); hi=max(ax,bx); flo=f(lo); fhi=f(hi)
    if(flo==0.0_dp)then;root=lo;return;end if
    if(fhi==0.0_dp)then;root=hi;return;end if
    if(flo*fhi>0.0_dp)then
      root=nan_dp(); if(present(status))status=ecd_invalid; return
    end if
    do it=1,nmax
      x_tol=t+2.0_dp*epsilon(1.0_dp)*max(abs(lo),abs(hi))
      if(abs(hi-lo)<=x_tol)then;root=0.5_dp*(lo+hi);return;end if
      ! Secant interpolation, safeguarded by bisection.  Reject a secant
      ! candidate that lies too close to an endpoint, where rounding can
      ! stall progress for very flat distribution tails.
      if(abs(fhi-flo)>tiny(1.0_dp))then
        c=(lo*fhi-hi*flo)/(fhi-flo)
      else
        c=0.5_dp*(lo+hi)
      end if
      if(c<=lo+0.1_dp*(hi-lo) .or. c>=hi-0.1_dp*(hi-lo))c=0.5_dp*(lo+hi)
      fc=f(c)
      if(fc==0.0_dp)then;root=c;return;end if
      if(flo*fc<0.0_dp)then
        hi=c; fhi=fc
      else
        lo=c; flo=fc
      end if
    end do
    root=0.5_dp*(lo+hi)
    if(present(status))status=ecd_no_convergence
  end function brent_root

  function bessel_k(nu, x) result(v)
    real(dp), intent(in) :: nu,x
    real(dp) :: v,k0,k1,k2,nur
    integer :: m,j
    if(x<0.0_dp) then; v=nan_dp(); return; end if
    if(x==0.0_dp) then; v=huge(1.0_dp); return; end if
    m=nint(nu+0.5_dp)
    if(abs(nu-(real(m,dp)-0.5_dp))<1.0e-12_dp .and. m>=1) then
      k0=sqrt(pi/(2.0_dp*x))*exp(-x)
      if(m==1) then; v=k0; return; end if
      k1=k0*(1.0_dp+1.0_dp/x)
      if(m==2) then; v=k1; return; end if
      do j=2,m-1
        nur=real(j,dp)-0.5_dp
        k2=k0+2.0_dp*nur/x*k1
        k0=k1; k1=k2
      end do
      v=k1
    else
      v=integrate_adaptive(kint,0.0_dp,20.0_dp,1.0e-9_dp,1.0e-12_dp)
    end if
  contains
    function kint(t) result(y)
      real(dp), intent(in) :: t
      real(dp) :: y,lg
      lg=-x*cosh(t)+log(max(cosh(nu*t),tiny(1.0_dp)))
      if(lg < log(tiny(1.0_dp))) then; y=0.0_dp; else; y=exp(lg); end if
    end function kint
  end function bessel_k

  function dawson_f(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y,ax,term,sumv
    integer :: n
    ax=abs(x)
    if(ax<0.2_dp) then
      term=x; sumv=x
      do n=1,40
        term=term*(-2.0_dp*x*x)/real(2*n+1,dp)
        sumv=sumv+term
        if(abs(term)<1e-16_dp*max(1.0_dp,abs(sumv))) exit
      end do
      y=sumv
    else if(ax>6.0_dp) then
      term=1.0_dp/(2.0_dp*x); sumv=term
      do n=1,12
        term=term*real(2*n-1,dp)/(2.0_dp*x*x)
        sumv=sumv+term
      end do
      y=sumv
    else
      y=exp(-x*x)*integrate_adaptive(exp_t2,0.0_dp,x,1e-11_dp,1e-13_dp)
    end if
  contains
    function exp_t2(t) result(z)
      real(dp), intent(in) :: t
      real(dp) :: z
      z=exp(t*t)
    end function exp_t2
  end function dawson_f

  function erfi_f(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    y=2.0_dp/sqrt_pi*exp(x*x)*dawson_f(x)
  end function erfi_f

  function erfcx_f(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y
    if(x<25.0_dp) then
      y=exp(x*x)*erfc(x)
    else
      y=(1.0_dp/x/sqrt_pi)*(1.0_dp-0.5_dp/x**2+0.75_dp/x**4-1.875_dp/x**6)
    end if
  end function erfcx_f

  function erfq(x,sgn) result(y)
    real(dp), intent(in) :: x
    integer, intent(in) :: sgn
    real(dp) :: y
    if(sgn==-1) then
      y=2.0_dp*dawson_f(x)
    else if(sgn==1) then
      y=sqrt_pi*erfcx_f(x)
    else
      y=nan_dp()
    end if
  end function erfq

  function erfq_sum(x,sgn) result(y)
    real(dp), intent(in) :: x
    integer, intent(in) :: sgn
    real(dp) :: y,term
    integer :: m,mx
    if((sgn/=1 .and. sgn/=-1) .or. x==0.0_dp) then; y=nan_dp(); return; end if
    mx=min(25,max(0,int(floor(x*x))))
    y=0.0_dp
    do m=0,mx
      term=gamma(real(m,dp)+0.5_dp)*real(-sgn,dp)**m/(x*x)**m
      y=y+term
    end do
    y=y/(sqrt_pi*x)
  end function erfq_sum

  function log_gamma_sign(x,sgn) result(lg)
    real(dp), intent(in) :: x
    integer, intent(out) :: sgn
    real(dp) :: lg
    if(x>0.0_dp) then
      lg=log_gamma(x); sgn=1
    else
      lg=log(pi)-log(abs(sin(pi*x)))-log_gamma(1.0_dp-x)
      sgn=merge(1,-1,sin(pi*x)>0.0_dp)
    end if
  end function log_gamma_sign

  function rational_approx(x,max_den) result(frac)
    real(dp), intent(in) :: x
    integer, intent(in), optional :: max_den
    integer :: frac(2), md,q,p
    real(dp) :: best,err
    md=10000; if(present(max_den)) md=max_den
    best=huge(1.0_dp); frac=[0,1]
    do q=1,md
      p=nint(x*real(q,dp)); err=abs(x-real(p,dp)/real(q,dp))
      if(err<best) then; best=err; frac=[p,q]; end if
      if(err<1e-12_dp) exit
    end do
  end function rational_approx

  function hypergeom_2f0(s,x,order) result(q)
    real(dp), intent(in) :: s,x
    integer, intent(in) :: order
    real(dp) :: q,si,xi
    integer :: i
    q=1.0_dp; si=1.0_dp; xi=1.0_dp
    do i=1,order
      si=si*(s-real(i,dp)); xi=xi*x; q=q+si/xi
    end do
  end function hypergeom_2f0

end module ecd_math
