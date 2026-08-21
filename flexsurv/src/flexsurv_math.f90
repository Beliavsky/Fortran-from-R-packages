! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_math
  use flexsurv_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, &
    ieee_quiet_nan, ieee_is_finite
  implicit none
  private

  real(dp), parameter, public :: fs_pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: sqrt2 = 1.41421356237309504880168872420969807857_dp
  real(dp), parameter :: log_sqrt_2pi = 0.918938533204672741780329736405617639861_dp

  public :: normal_pdf, normal_logpdf, normal_cdf, normal_logcdf, normal_quantile
  public :: regularized_gamma_p, regularized_gamma_q, gamma_quantile
  public :: regularized_beta, beta_quantile
  public :: log_beta, log1mexp, logdiffexp, safe_log, log1p_fs, expm1_fs
  public :: rng_seed, rng_normal, rng_gamma, rng_beta
  public :: integrate_simpson, integrate_gauss_legendre
  public :: brent_root, near_positive_definite
  public :: logsumexp

  abstract interface
    function scalar_fn(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fn
  end interface

contains

  pure real(dp) function log1p_fs(x) result(y)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-8_dp) then
      y = x - 0.5_dp*x*x + x**3/3.0_dp - 0.25_dp*x**4 + 0.2_dp*x**5
    else
      y = log(1.0_dp + x)
    end if
  end function log1p_fs

  pure real(dp) function expm1_fs(x) result(y)
    real(dp), intent(in) :: x
    if (abs(x) < 1.0e-8_dp) then
      y = x + 0.5_dp*x*x + x**3/6.0_dp + x**4/24.0_dp + x**5/120.0_dp
    else
      y = exp(x) - 1.0_dp
    end if
  end function expm1_fs

  pure real(dp) function safe_log(x) result(y)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      y = log(x)
    else if (x == 0.0_dp) then
      y = -ieee_value(0.0_dp, ieee_positive_inf)
    else
      y = ieee_value(0.0_dp, ieee_quiet_nan)
    end if
  end function safe_log

  pure real(dp) function log1mexp(a) result(y)
    real(dp), intent(in) :: a
    if (a > 0.0_dp) then
      y = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (a < -0.693147180559945309417232121458176568_dp) then
      y = log1p_fs(-exp(a))
    else if (a < 0.0_dp) then
      y = log(-expm1_fs(a))
    else
      y = -ieee_value(0.0_dp, ieee_positive_inf)
    end if
  end function log1mexp

  pure real(dp) function logdiffexp(a, b) result(y)
    real(dp), intent(in) :: a, b
    if (b > a) then
      y = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (.not. ieee_is_finite(b)) then
      y = a
    else
      y = a + log1mexp(b - a)
    end if
  end function logdiffexp

  pure real(dp) function logsumexp(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) == 0) then
      y = -ieee_value(0.0_dp, ieee_positive_inf)
      return
    end if
    m = maxval(x)
    if (.not. ieee_is_finite(m)) then
      y = m
    else
      y = m + log(sum(exp(x - m)))
    end if
  end function logsumexp

  pure real(dp) function normal_logpdf(x) result(y)
    real(dp), intent(in) :: x
    y = -0.5_dp*x*x - log_sqrt_2pi
  end function normal_logpdf

  pure real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(normal_logpdf(x))
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    y = 0.5_dp * erfc(-x/sqrt2)
  end function normal_cdf

  pure real(dp) function normal_logcdf(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: p
    if (x < -10.0_dp) then
      ! Mills expansion, enough for survival likelihood tails.
      y = normal_logpdf(x) - log(-x) + log(1.0_dp - 1.0_dp/(x*x) + &
        3.0_dp/(x**4) - 15.0_dp/(x**6))
    else
      p = normal_cdf(x)
      y = log(p)
    end if
  end function normal_logcdf

  pure real(dp) function normal_quantile(p) result(x)
    ! Acklam rational approximation.
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
      -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
      -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
      -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
      -1.328068155288572e+01_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
      -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
       4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
       2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
    real(dp) :: q, r, e, u
    if (p <= 0.0_dp) then
      x = -ieee_value(0.0_dp, ieee_positive_inf)
      return
    else if (p >= 1.0_dp) then
      x = ieee_value(0.0_dp, ieee_positive_inf)
      return
    end if
    if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
    ! Halley refinement.
    e = normal_cdf(x) - p
    u = e / normal_pdf(x)
    x = x - u/(1.0_dp + 0.5_dp*x*u)
  end function normal_quantile

  pure real(dp) function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: itmax = 500
    real(dp), parameter :: eps = 2.0e-15_dp, fpmin = 1.0e-300_dp
    real(dp) :: sumv, del, ap, b, c, d, h, an, gln
    integer :: n
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x == 0.0_dp) then
      p = 0.0_dp
      return
    end if
    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, itmax
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x + a*log(x) - gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n = 1, itmax
        an = -real(n,dp)*(real(n,dp)-a)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del - 1.0_dp) <= eps) exit
      end do
      p = 1.0_dp - exp(-x + a*log(x) - gln)*h
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function regularized_gamma_p

  pure real(dp) function regularized_gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    q = 1.0_dp - regularized_gamma_p(a, x)
    q = max(0.0_dp, min(1.0_dp, q))
  end function regularized_gamma_q

  pure real(dp) function log_beta(a, b) result(y)
    real(dp), intent(in) :: a, b
    y = log_gamma(a) + log_gamma(b) - log_gamma(a+b)
  end function log_beta

  pure real(dp) function betacf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 500
    real(dp), parameter :: eps=3.0e-15_dp, fpmin=1.0e-300_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp
    d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin)d=fpmin
    d=1.0_dp/d; h=d
    do m=1,maxit
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/ &
         ((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps)exit
    end do
    cf=h
  end function betacf

  pure real(dp) function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x == 0.0_dp) then
      p = 0.0_dp
    else if (x == 1.0_dp) then
      p = 1.0_dp
    else
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log1p_fs(-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        p = bt*betacf(a,b,x)/a
      else
        p = 1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
      end if
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function regularized_beta

  real(dp) function gamma_quantile(p, shape, scale) result(x)
    real(dp), intent(in) :: p, shape
    real(dp), intent(in), optional :: scale
    real(dp) :: lo, hi, mid, sc
    integer :: it
    sc=1.0_dp; if(present(scale))sc=scale
    if(p<=0.0_dp)then; x=0.0_dp; return; end if
    if(p>=1.0_dp)then; x=ieee_value(0.0_dp,ieee_positive_inf); return; end if
    lo=0.0_dp
    hi=max(1.0_dp, shape + 10.0_dp*sqrt(max(shape,1.0_dp)))
    do while(regularized_gamma_p(shape,hi)<p)
      hi=2.0_dp*hi
      if(hi>1.0e300_dp)exit
    end do
    do it=1,180
      mid=0.5_dp*(lo+hi)
      if(regularized_gamma_p(shape,mid)<p)then; lo=mid; else; hi=mid; end if
    end do
    x=0.5_dp*(lo+hi)*sc
  end function gamma_quantile

  real(dp) function beta_quantile(p, a, b) result(x)
    real(dp), intent(in) :: p,a,b
    real(dp)::lo,hi,mid
    integer::it
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=1.0_dp;return;end if
    lo=0.0_dp;hi=1.0_dp
    do it=1,160
      mid=0.5_dp*(lo+hi)
      if(regularized_beta(mid,a,b)<p)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function beta_quantile

  subroutine rng_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=mod(abs(seed)+104729*i+37*i*i, huge(1)-1)+1
    end do
    call random_seed(put=put)
  end subroutine rng_seed

  real(dp) function rng_normal() result(z)
    real(dp) :: u1,u2
    call random_number(u1); call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*fs_pi*u2)
  end function rng_normal

  recursive real(dp) function rng_gamma(shape, scale) result(x)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::scale
    real(dp)::d,c,z,u,sc
    sc=1.0_dp;if(present(scale))sc=scale
    if(shape<=0.0_dp)then;x=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(shape<1.0_dp)then
      call random_number(u)
      x=rng_gamma(shape+1.0_dp)*u**(1.0_dp/shape)*sc
      return
    end if
    d=shape-1.0_dp/3.0_dp;c=1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z=rng_normal()
        if(1.0_dp+c*z>0.0_dp)exit
      end do
      x=(1.0_dp+c*z)**3
      call random_number(u)
      if(u<1.0_dp-0.0331_dp*z**4)exit
      if(log(u)<0.5_dp*z*z+d*(1.0_dp-x+log(x)))exit
    end do
    x=d*x*sc
  end function rng_gamma

  real(dp) function rng_beta(a,b) result(x)
    real(dp),intent(in)::a,b
    real(dp)::g1,g2
    g1=rng_gamma(a);g2=rng_gamma(b);x=g1/(g1+g2)
  end function rng_beta

  real(dp) function integrate_simpson(fn,a,b,tol,max_depth) result(val)
    procedure(scalar_fn)::fn
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::max_depth
    real(dp)::eps,fa,fb,fc,s
    integer::md
    eps=1.0e-9_dp;if(present(tol))eps=tol
    md=20;if(present(max_depth))md=max_depth
    if(a==b)then;val=0.0_dp;return;end if
    fa=fn(a);fb=fn(b);fc=fn(0.5_dp*(a+b))
    s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    val=adaptive_simpson(fn,a,b,fa,fb,fc,s,eps,md)
  end function integrate_simpson

  recursive real(dp) function adaptive_simpson(fn,a,b,fa,fb,fc,s,eps,depth) result(v)
    procedure(scalar_fn)::fn
    real(dp),intent(in)::a,b,fa,fb,fc,s,eps
    integer,intent(in)::depth
    real(dp)::c,d,e,fd,fe,s1,s2
    c=0.5_dp*(a+b);d=0.5_dp*(a+c);e=0.5_dp*(c+b)
    fd=fn(d);fe=fn(e)
    s1=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
    s2=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
    if(depth<=0 .or. abs(s1+s2-s)<=15.0_dp*eps)then
      v=s1+s2+(s1+s2-s)/15.0_dp
    else
      v=adaptive_simpson(fn,a,c,fa,fc,fd,s1,0.5_dp*eps,depth-1)+ &
        adaptive_simpson(fn,c,b,fc,fb,fe,s2,0.5_dp*eps,depth-1)
    end if
  end function adaptive_simpson

  real(dp) function integrate_gauss_legendre(fn,a,b,n) result(val)
    procedure(scalar_fn)::fn
    real(dp),intent(in)::a,b
    integer,intent(in),optional::n
    integer::nn,i,j,m
    real(dp)::z,z1,p1,p2,p3,pp,xm,xl,w
    nn=64;if(present(n))nn=n
    m=(nn+1)/2;xm=0.5_dp*(b+a);xl=0.5_dp*(b-a);val=0.0_dp
    do i=1,m
      z=cos(fs_pi*(real(i,dp)-0.25_dp)/(real(nn,dp)+0.5_dp))
      do
        p1=1.0_dp;p2=0.0_dp
        do j=1,nn
          p3=p2;p2=p1
          p1=((2.0_dp*real(j,dp)-1.0_dp)*z*p2-(real(j,dp)-1.0_dp)*p3)/real(j,dp)
        end do
        pp=real(nn,dp)*(z*p1-p2)/(z*z-1.0_dp)
        z1=z;z=z1-p1/pp
        if(abs(z-z1)<2.0e-15_dp)exit
      end do
      w=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
      val=val+w*(fn(xm-xl*z)+fn(xm+xl*z))
    end do
  end function integrate_gauss_legendre

  real(dp) function brent_root(fn,a,b,tol,status) result(root)
    procedure(scalar_fn)::fn
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::tol
    integer,intent(out),optional::status
    real(dp)::lo,hi,mid,flo,fhi,fmid,eps
    integer::it
    eps=1.0e-10_dp;if(present(tol))eps=tol
    lo=a;hi=b;flo=fn(lo);fhi=fn(hi)
    if(flo==0.0_dp)then;root=lo;if(present(status))status=0;return;end if
    if(fhi==0.0_dp)then;root=hi;if(present(status))status=0;return;end if
    if(flo*fhi>0.0_dp)then
      root=ieee_value(0.0_dp,ieee_quiet_nan);if(present(status))status=1;return
    end if
    do it=1,300
      mid=0.5_dp*(lo+hi);fmid=fn(mid)
      if(abs(fmid)<eps .or. abs(hi-lo)<eps*max(1.0_dp,abs(mid)))exit
      if(flo*fmid<=0.0_dp)then;hi=mid;fhi=fmid;else;lo=mid;flo=fmid;end if
    end do
    root=0.5_dp*(lo+hi);if(present(status))status=0
  end function brent_root

  subroutine near_positive_definite(a, out, eps)
    ! Self-contained Higham-style diagonal-loading repair.  It preserves
    ! symmetry and is sufficient for covariance simulation in this port.
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::out(size(a,1),size(a,2))
    real(dp),intent(in),optional::eps
    real(dp)::jitter,base
    integer::n,i,k
    real(dp),allocatable::l(:,:)
    n=size(a,1);out=0.5_dp*(a+transpose(a));allocate(l(n,n))
    base=1.0e-10_dp;if(present(eps))base=eps
    jitter=0.0_dp
    do k=0,14
      l=out
      do i=1,n;l(i,i)=l(i,i)+jitter;end do
      if(chol_ok(l))then
        if(jitter>0.0_dp)then
          do i=1,n;out(i,i)=out(i,i)+jitter;end do
        end if
        return
      end if
      if(k==0)then;jitter=base*max(1.0_dp,maxval(abs(out)));else;jitter=10.0_dp*jitter;end if
    end do
  end subroutine near_positive_definite

  logical function chol_ok(a) result(ok)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable::l(:,:)
    real(dp)::s
    integer::n,i,j,k
    n=size(a,1);allocate(l(n,n));l=0.0_dp;ok=.true.
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1;s=s-l(i,k)*l(j,k);end do
        if(i==j)then
          if(s<=0.0_dp .or. .not.ieee_is_finite(s))then;ok=.false.;return;end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end function chol_ok

end module flexsurv_math
