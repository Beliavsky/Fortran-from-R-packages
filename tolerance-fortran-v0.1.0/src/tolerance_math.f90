! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_math
  use tolerance_kinds, only : dp, pi, sqrt2, huge_dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private

  public :: normal_pdf, normal_cdf, normal_quantile
  public :: gamma_p, gamma_q, beta_i, beta_quantile, digamma, polygamma
  public :: gamma_quantile, chisq_cdf, chisq_quantile
  public :: student_t_cdf, student_t_quantile, noncentral_t_cdf, noncentral_t_quantile
  public :: f_cdf, f_quantile, noncentral_chisq_cdf, noncentral_chisq_quantile
  public :: binom_pmf, binom_cdf, binom_quantile
  public :: poisson_pmf, poisson_cdf, poisson_quantile
  public :: negbin_pmf, negbin_cdf, negbin_quantile
  public :: hypergeom_pmf, hypergeom_cdf, hypergeom_quantile
  public :: log_choose, log_beta_fn, clamp
  public :: sample_mean, sample_variance, sample_sd, sample_median, sample_quantile
  public :: rng_normal, rng_gamma, rng_beta, rng_chisq, rng_poisson, rng_binomial
  public :: adaptive_simpson, bisect_root
  public :: solve_linear, invert_matrix, cholesky_lower, solve_spd, invert_spd
  public :: mvn_random

  abstract interface
    function scalar_function(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_function
  end interface

contains

  elemental real(dp) function clamp(x, lo, hi)
    real(dp), intent(in) :: x, lo, hi
    clamp = max(lo, min(hi, x))
  end function clamp

  elemental real(dp) function normal_pdf(x)
    real(dp), intent(in) :: x
    normal_pdf = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  elemental real(dp) function normal_cdf(x)
    real(dp), intent(in) :: x
    normal_cdf = 0.5_dp*erfc(-x/sqrt2)
  end function normal_cdf

  elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376e1_dp, a2=2.209460984245205e2_dp
    real(dp), parameter :: a3=-2.759285104469687e2_dp, a4=1.383577518672690e2_dp
    real(dp), parameter :: a5=-3.066479806614716e1_dp, a6=2.506628277459239_dp
    real(dp), parameter :: b1=-5.447609879822406e1_dp, b2=1.615858368580409e2_dp
    real(dp), parameter :: b3=-1.556989798598866e2_dp, b4=6.680131188771972e1_dp
    real(dp), parameter :: b5=-1.328068155288572e1_dp
    real(dp), parameter :: c1=-7.784894002430293e-3_dp, c2=-3.223964580411365e-1_dp
    real(dp), parameter :: c3=-2.400758277161838_dp, c4=-2.549732539343734_dp
    real(dp), parameter :: c5=4.374664141464968_dp, c6=2.938163982698783_dp
    real(dp), parameter :: d1=7.784695709041462e-3_dp, d2=3.224671290700398e-1_dp
    real(dp), parameter :: d3=2.445134137142996_dp, d4=3.754408661907416_dp
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r, e
    integer :: iter

    if (p <= 0.0_dp) then
      x = -huge_dp
      return
    else if (p >= 1.0_dp) then
      x = huge_dp
      return
    end if
    if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
          ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
          (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
           ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
    ! Two Halley corrections improve the Acklam approximation to near machine precision.
    do iter=1,2
      e = normal_cdf(x)-p
      x = x - e/(normal_pdf(x) + 0.5_dp*x*e)
    end do
  end function normal_quantile

  elemental real(dp) function digamma(x) result(v)
    real(dp),intent(in)::x
    real(dp)::y,r
    if(x<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    y=x;r=0.0_dp
    do while(y<10.0_dp)
      r=r-1.0_dp/y;y=y+1.0_dp
    end do
    v=r+log(y)-0.5_dp/y-1.0_dp/(12.0_dp*y*y)+1.0_dp/(120.0_dp*y**4) - &
      1.0_dp/(252.0_dp*y**6)+1.0_dp/(240.0_dp*y**8)-5.0_dp/(660.0_dp*y**10)
  end function digamma

  real(dp) function polygamma(order,x) result(v)
    integer,intent(in)::order
    real(dp),intent(in)::x
    real(dp)::y,r,sgn,fact,z
    integer::j
    if(order<0 .or. x<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    if(order==0)then;v=digamma(x);return;end if
    y=x;r=0.0_dp;fact=1.0_dp
    do j=2,order;fact=fact*real(j,dp);end do
    sgn=merge(1.0_dp,-1.0_dp,mod(order+1,2)==0)
    do while(y<12.0_dp)
      ! psi_m(y) = psi_m(y+1) - (-1)^m m!/y^(m+1)
      r=r-merge(1.0_dp,-1.0_dp,mod(order,2)==0)*fact/y**(order+1)
      y=y+1.0_dp
    end do
    z=hurwitz_zeta_int(order+1,y)
    v=r+sgn*fact*z
  end function polygamma

  real(dp) function hurwitz_zeta_int(s,a) result(z)
    integer,intent(in)::s
    real(dp),intent(in)::a
    integer,parameter::n=24
    integer::k
    real(dp)::x,ss
    ss=real(s,dp);z=0.0_dp
    do k=0,n-1;z=z+(a+real(k,dp))**(-ss);end do
    x=a+real(n,dp)
    z=z+x**(1.0_dp-ss)/(ss-1.0_dp)+0.5_dp*x**(-ss)+ss*x**(-ss-1.0_dp)/12.0_dp - &
      ss*(ss+1.0_dp)*(ss+2.0_dp)*x**(-ss-3.0_dp)/720.0_dp + &
      ss*(ss+1.0_dp)*(ss+2.0_dp)*(ss+3.0_dp)*(ss+4.0_dp)*x**(-ss-5.0_dp)/30240.0_dp
  end function hurwitz_zeta_int

  pure real(dp) function gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: itmax=10000
    real(dp), parameter :: eps=2.0e-15_dp, fpmin=tiny(1.0_dp)/eps
    integer :: n
    real(dp) :: ap, del, sumv, b, c, d, h, an, q
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x == 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x < a+1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n=1,itmax
        ap = ap+1.0_dp
        del = del*x/ap
        sumv = sumv+del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x+1.0_dp-a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n=1,itmax
        an = -real(n,dp)*(real(n,dp)-a)
        b = b+2.0_dp
        d = an*d+b
        if (abs(d) < fpmin) d=fpmin
        c = b+an/c
        if (abs(c) < fpmin) c=fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      q = exp(-x+a*log(x)-log_gamma(a))*h
      p = 1.0_dp-q
    end if
    p = clamp(p,0.0_dp,1.0_dp)
  end function gamma_p

  pure real(dp) function gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    q = 1.0_dp-gamma_p(a,x)
    q = clamp(q,0.0_dp,1.0_dp)
  end function gamma_q

  pure real(dp) function beta_i(x, a, b) result(v)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      v = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      v=0.0_dp; return
    else if (x >= 1.0_dp) then
      v=1.0_dp; return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      v = bt*beta_cf(a,b,x)/a
    else
      v = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    v=clamp(v,0.0_dp,1.0_dp)
  end function beta_i

  pure real(dp) function beta_cf(a,b,x) result(h)
    real(dp), intent(in) :: a,b,x
    integer, parameter :: maxit=10000
    real(dp), parameter :: eps=2.0e-15_dp, fpmin=tiny(1.0_dp)/eps
    integer :: m, m2
    real(dp) :: qab,qap,qam,c,d,aa,del
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp
    d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin)d=fpmin
    d=1.0_dp/d
    h=d
    do m=1,maxit
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
      c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
      d=1.0_dp/d; del=d*c; h=h*del
      if(abs(del-1.0_dp)<=eps)exit
    end do
  end function beta_cf

  elemental real(dp) function log_beta_fn(a,b)
    real(dp), intent(in) :: a,b
    log_beta_fn=log_gamma(a)+log_gamma(b)-log_gamma(a+b)
  end function log_beta_fn

  elemental real(dp) function log_choose(n,k)
    integer, intent(in) :: n,k
    if(k<0 .or. k>n)then
      log_choose=-huge_dp
    else
      log_choose=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
    end if
  end function log_choose

  real(dp) function beta_quantile(p,a,b) result(x)
    real(dp),intent(in)::p,a,b
    real(dp)::lo,hi,mid
    integer::i
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=1.0_dp;return;end if
    lo=0.0_dp;hi=1.0_dp
    do i=1,160
      mid=0.5_dp*(lo+hi)
      if(beta_i(mid,a,b)<p)then;lo=mid;else;hi=mid;end if
      if(hi-lo<2.0e-14_dp)exit
    end do
    x=0.5_dp*(lo+hi)
  end function beta_quantile

  real(dp) function gamma_quantile(p,a,scale) result(x)
    real(dp), intent(in)::p,a
    real(dp), intent(in), optional::scale
    real(dp)::lo,hi,s
    integer::i
    s=1.0_dp; if(present(scale))s=scale
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=huge_dp;return;end if
    lo=0.0_dp; hi=max(1.0_dp,a*s)
    do while(gamma_p(a,hi/s)<p .and. hi<huge_dp/4.0_dp)
      hi=2.0_dp*hi
    end do
    do i=1,160
      x=0.5_dp*(lo+hi)
      if(gamma_p(a,x/s)<p)then;lo=x;else;hi=x;end if
      if(abs(hi-lo)<=1.0e-12_dp*max(1.0_dp,x))exit
    end do
    x=0.5_dp*(lo+hi)
  end function gamma_quantile

  elemental real(dp) function chisq_cdf(x,df)
    real(dp),intent(in)::x,df
    if(x<=0.0_dp)then
      chisq_cdf=0.0_dp
    else
      chisq_cdf=gamma_p(0.5_dp*df,0.5_dp*x)
    end if
  end function chisq_cdf

  real(dp) function chisq_quantile(p,df) result(x)
    real(dp),intent(in)::p,df
    x=gamma_quantile(p,0.5_dp*df,2.0_dp)
  end function chisq_quantile

  elemental real(dp) function student_t_cdf(x,df) result(p)
    real(dp),intent(in)::x,df
    real(dp)::z
    if(df<=0.0_dp)then;p=ieee_value(0.0_dp, ieee_quiet_nan);return;end if
    if(x==0.0_dp)then;p=0.5_dp;return;end if
    z=df/(df+x*x)
    if(x>0.0_dp)then
      p=1.0_dp-0.5_dp*beta_i(z,0.5_dp*df,0.5_dp)
    else
      p=0.5_dp*beta_i(z,0.5_dp*df,0.5_dp)
    end if
  end function student_t_cdf

  real(dp) function student_t_quantile(p,df) result(x)
    real(dp),intent(in)::p,df
    real(dp)::lo,hi,mid
    integer::i
    if(p<=0.0_dp)then;x=-huge_dp;return;end if
    if(p>=1.0_dp)then;x=huge_dp;return;end if
    lo=-1.0_dp;hi=1.0_dp
    do while(student_t_cdf(lo,df)>p);lo=2.0_dp*lo;end do
    do while(student_t_cdf(hi,df)<p);hi=2.0_dp*hi;end do
    do i=1,180
      mid=0.5_dp*(lo+hi)
      if(student_t_cdf(mid,df)<p)then;lo=mid;else;hi=mid;end if
      if(abs(hi-lo)<1.0e-12_dp*max(1.0_dp,abs(mid)))exit
    end do
    x=0.5_dp*(lo+hi)
  end function student_t_quantile

  elemental real(dp) function f_cdf(x,df1,df2) result(p)
    real(dp),intent(in)::x,df1,df2
    if(x<=0.0_dp)then
      p=0.0_dp
    else
      p=beta_i(df1*x/(df1*x+df2),0.5_dp*df1,0.5_dp*df2)
    end if
  end function f_cdf

  real(dp) function f_quantile(p,df1,df2) result(x)
    real(dp),intent(in)::p,df1,df2
    real(dp)::lo,hi,mid
    integer::i
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=huge_dp;return;end if
    lo=0.0_dp;hi=1.0_dp
    do while(f_cdf(hi,df1,df2)<p);hi=2.0_dp*hi;end do
    do i=1,180
      mid=0.5_dp*(lo+hi)
      if(f_cdf(mid,df1,df2)<p)then;lo=mid;else;hi=mid;end if
      if(abs(hi-lo)<1.0e-12_dp*max(1.0_dp,mid))exit
    end do
    x=0.5_dp*(lo+hi)
  end function f_quantile

  real(dp) function noncentral_chisq_cdf(x,df,ncp) result(p)
    real(dp),intent(in)::x,df,ncp
    real(dp)::lambda,w,sumw,term
    integer::k
    if(x<=0.0_dp)then;p=0.0_dp;return;end if
    if(ncp<=1.0e-14_dp)then;p=chisq_cdf(x,df);return;end if
    lambda=0.5_dp*ncp
    if(lambda<350.0_dp)then
      w=exp(-lambda);sumw=w
      p=w*chisq_cdf(x,df)
      do k=1,100000
        w=w*lambda/real(k,dp)
        term=w*chisq_cdf(x,df+2.0_dp*real(k,dp))
        p=p+term;sumw=sumw+w
        if(w<1.0e-15_dp*max(sumw,1.0_dp) .and. real(k,dp)>lambda+10.0_dp*sqrt(lambda))exit
      end do
    else
      ! Normal approximation is used only for extremely large noncentrality.
      p=normal_cdf((x-(df+ncp))/sqrt(2.0_dp*(df+2.0_dp*ncp)))
    end if
    p=clamp(p,0.0_dp,1.0_dp)
  end function noncentral_chisq_cdf

  real(dp) function noncentral_chisq_quantile(p,df,ncp) result(x)
    real(dp),intent(in)::p,df,ncp
    real(dp)::lo,hi,mid
    integer::i
    if(p<=0.0_dp)then;x=0.0_dp;return;end if
    if(p>=1.0_dp)then;x=huge_dp;return;end if
    lo=0.0_dp;hi=max(1.0_dp,df+ncp+10.0_dp*sqrt(2.0_dp*(df+2.0_dp*ncp)))
    do while(noncentral_chisq_cdf(hi,df,ncp)<p);hi=2.0_dp*hi;end do
    do i=1,180
      mid=0.5_dp*(lo+hi)
      if(noncentral_chisq_cdf(mid,df,ncp)<p)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function noncentral_chisq_quantile

  real(dp) function noncentral_t_cdf(t,df,ncp) result(p)
    real(dp),intent(in)::t,df,ncp
    real(dp)::upper,a
    if(df<=0.0_dp)then;p=ieee_value(0.0_dp, ieee_quiet_nan);return;end if
    if(abs(ncp)<1.0e-14_dp)then;p=student_t_cdf(t,df);return;end if
    a=0.5_dp*df
    upper=gamma_quantile(1.0_dp-1.0e-12_dp,a,1.0_dp)
    p=adaptive_simpson(integrand,0.0_dp,upper,1.0e-10_dp,24)
    p=clamp(p,0.0_dp,1.0_dp)
  contains
    function integrand(y) result(v)
      real(dp),intent(in)::y
      real(dp)::v,lg
      if(y<=0.0_dp)then
        if(a<1.0_dp)then
          v=0.0_dp
        else if(a==1.0_dp)then
          v=normal_cdf(-ncp)
        else
          v=0.0_dp
        end if
      else
        lg=(a-1.0_dp)*log(y)-y-log_gamma(a)
        v=normal_cdf(t*sqrt(2.0_dp*y/df)-ncp)*exp(lg)
      end if
    end function integrand
  end function noncentral_t_cdf

  real(dp) function noncentral_t_quantile(p,df,ncp) result(x)
    real(dp),intent(in)::p,df,ncp
    real(dp)::lo,hi,mid
    integer::i
    if(p<=0.0_dp)then;x=-huge_dp;return;end if
    if(p>=1.0_dp)then;x=huge_dp;return;end if
    lo=min(-1.0_dp,ncp-10.0_dp);hi=max(1.0_dp,ncp+10.0_dp)
    do while(noncentral_t_cdf(lo,df,ncp)>p);lo=2.0_dp*lo-1.0_dp;end do
    do while(noncentral_t_cdf(hi,df,ncp)<p);hi=2.0_dp*hi+1.0_dp;end do
    do i=1,100
      mid=0.5_dp*(lo+hi)
      if(noncentral_t_cdf(mid,df,ncp)<p)then;lo=mid;else;hi=mid;end if
      if(abs(hi-lo)<5.0e-10_dp*max(1.0_dp,abs(mid)))exit
    end do
    x=0.5_dp*(lo+hi)
  end function noncentral_t_quantile

  elemental real(dp) function binom_pmf(k,n,p) result(v)
    integer,intent(in)::k,n
    real(dp),intent(in)::p
    if(k<0 .or. k>n .or. p<0.0_dp .or. p>1.0_dp)then;v=0.0_dp;return;end if
    if(p==0.0_dp)then;v=merge(1.0_dp,0.0_dp,k==0);return;end if
    if(p==1.0_dp)then;v=merge(1.0_dp,0.0_dp,k==n);return;end if
    v=exp(log_choose(n,k)+real(k,dp)*log(p)+real(n-k,dp)*log(1.0_dp-p))
  end function binom_pmf

  elemental real(dp) function binom_cdf(k,n,p) result(v)
    integer,intent(in)::k,n
    real(dp),intent(in)::p
    if(k<0)then;v=0.0_dp;return;end if
    if(k>=n)then;v=1.0_dp;return;end if
    v=beta_i(1.0_dp-p,real(n-k,dp),real(k+1,dp))
  end function binom_cdf

  integer function binom_quantile(prob,n,p) result(k)
    real(dp),intent(in)::prob,p
    integer,intent(in)::n
    integer::lo,hi,mid
    if(prob<=0.0_dp)then;k=0;return;end if
    if(prob>=1.0_dp)then;k=n;return;end if
    lo=0;hi=n
    do while(lo<hi)
      mid=(lo+hi)/2
      if(binom_cdf(mid,n,p)>=prob)then;hi=mid;else;lo=mid+1;end if
    end do
    k=lo
  end function binom_quantile

  elemental real(dp) function poisson_pmf(k,lambda) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::lambda
    if(k<0 .or. lambda<0.0_dp)then;v=0.0_dp;return;end if
    if(lambda==0.0_dp)then;v=merge(1.0_dp,0.0_dp,k==0);return;end if
    v=exp(-lambda+real(k,dp)*log(lambda)-log_gamma(real(k+1,dp)))
  end function poisson_pmf

  elemental real(dp) function poisson_cdf(k,lambda) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::lambda
    if(k<0)then;v=0.0_dp;else;v=gamma_q(real(k+1,dp),lambda);end if
  end function poisson_cdf

  integer function poisson_quantile(prob,lambda) result(k)
    real(dp),intent(in)::prob,lambda
    integer::lo,hi,mid
    if(prob<=0.0_dp .or. lambda<=0.0_dp)then;k=0;return;end if
    lo=0;hi=max(1,int(lambda+10.0_dp*sqrt(lambda+1.0_dp)+20.0_dp))
    do while(poisson_cdf(hi,lambda)<prob .and. hi<huge(hi)/4);hi=2*hi;end do
    do while(lo<hi)
      mid=(lo+hi)/2
      if(poisson_cdf(mid,lambda)>=prob)then;hi=mid;else;lo=mid+1;end if
    end do
    k=lo
  end function poisson_quantile

  elemental real(dp) function negbin_pmf(k,size,p) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::size,p
    if(k<0 .or. size<=0.0_dp .or. p<=0.0_dp .or. p>1.0_dp)then;v=0.0_dp;return;end if
    v=exp(log_gamma(real(k,dp)+size)-log_gamma(size)-log_gamma(real(k+1,dp)) + &
      size*log(p)+real(k,dp)*log(1.0_dp-p))
  end function negbin_pmf

  elemental real(dp) function negbin_cdf(k,size,p) result(v)
    integer,intent(in)::k
    real(dp),intent(in)::size,p
    if(k<0)then;v=0.0_dp;else;v=beta_i(p,size,real(k+1,dp));end if
  end function negbin_cdf

  integer function negbin_quantile(prob,size,p) result(k)
    real(dp),intent(in)::prob,size,p
    integer::lo,hi,mid
    real(dp)::mu,var
    if(prob<=0.0_dp)then;k=0;return;end if
    if(p<=0.0_dp)then;k=huge(k);return;end if
    mu=size*(1.0_dp-p)/p;var=size*(1.0_dp-p)/(p*p)
    lo=0;hi=max(1,int(mu+12.0_dp*sqrt(var+1.0_dp)+20.0_dp))
    do while(negbin_cdf(hi,size,p)<prob .and. hi<huge(hi)/4);hi=2*hi;end do
    do while(lo<hi)
      mid=(lo+hi)/2
      if(negbin_cdf(mid,size,p)>=prob)then;hi=mid;else;lo=mid+1;end if
    end do
    k=lo
  end function negbin_quantile

  elemental real(dp) function hypergeom_pmf(k,good,bad,draws) result(v)
    integer,intent(in)::k,good,bad,draws
    integer::lo,hi
    lo=max(0,draws-bad);hi=min(draws,good)
    if(k<lo .or. k>hi)then;v=0.0_dp;return;end if
    v=exp(log_choose(good,k)+log_choose(bad,draws-k)-log_choose(good+bad,draws))
  end function hypergeom_pmf

  real(dp) function hypergeom_cdf(k,good,bad,draws) result(v)
    integer,intent(in)::k,good,bad,draws
    integer::j,lo,hi
    lo=max(0,draws-bad);hi=min(k,min(draws,good));v=0.0_dp
    if(hi<lo)return
    do j=lo,hi;v=v+hypergeom_pmf(j,good,bad,draws);end do
    v=clamp(v,0.0_dp,1.0_dp)
  end function hypergeom_cdf

  integer function hypergeom_quantile(prob,good,bad,draws) result(k)
    real(dp),intent(in)::prob
    integer,intent(in)::good,bad,draws
    integer::lo,hi,mid
    lo=max(0,draws-bad);hi=min(draws,good)
    if(prob<=0.0_dp)then;k=lo;return;end if
    do while(lo<hi)
      mid=(lo+hi)/2
      if(hypergeom_cdf(mid,good,bad,draws)>=prob)then;hi=mid;else;lo=mid+1;end if
    end do
    k=lo
  end function hypergeom_quantile

  real(dp) function sample_mean(x) result(v)
    real(dp),intent(in)::x(:)
    if(size(x)==0)then;v=0.0_dp;else;v=sum(x)/real(size(x),dp);end if
  end function sample_mean

  real(dp) function sample_variance(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::m
    if(size(x)<2)then;v=0.0_dp;return;end if
    m=sample_mean(x);v=sum((x-m)**2)/real(size(x)-1,dp)
  end function sample_variance

  real(dp) function sample_sd(x) result(v)
    real(dp),intent(in)::x(:)
    v=sqrt(max(0.0_dp,sample_variance(x)))
  end function sample_sd

  real(dp) function sample_median(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp),allocatable::y(:)
    integer::n
    y=x;call sort_real(y);n=size(y)
    if(mod(n,2)==1)then;v=y((n+1)/2);else;v=0.5_dp*(y(n/2)+y(n/2+1));end if
  end function sample_median

  real(dp) function sample_quantile(x,p) result(v)
    real(dp),intent(in)::x(:),p
    real(dp),allocatable::y(:)
    real(dp)::h,frac
    integer::n,j
    y=x;call sort_real(y);n=size(y)
    if(n==0)then;v=ieee_value(0.0_dp, ieee_quiet_nan);return;end if
    if(p<=0.0_dp)then;v=y(1);return;end if
    if(p>=1.0_dp)then;v=y(n);return;end if
    h=1.0_dp+real(n-1,dp)*p;j=floor(h);frac=h-real(j,dp)
    if(j>=n)then;v=y(n);else;v=(1.0_dp-frac)*y(j)+frac*y(j+1);end if
  end function sample_quantile

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  real(dp) function rng_normal() result(z)
    real(dp)::u1,u2
    call random_number(u1);call random_number(u2)
    u1=max(u1,tiny(1.0_dp));z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rng_normal

  recursive real(dp) function rng_gamma(shape,scale) result(x)
    real(dp),intent(in)::shape
    real(dp),intent(in),optional::scale
    real(dp)::d,c,z,u,s
    s=1.0_dp;if(present(scale))s=scale
    if(shape<=0.0_dp)then;x=ieee_value(0.0_dp, ieee_quiet_nan);return;end if
    if(shape<1.0_dp)then
      call random_number(u)
      x=rng_gamma(shape+1.0_dp,1.0_dp)*u**(1.0_dp/shape)*s
      return
    end if
    d=shape-1.0_dp/3.0_dp;c=1.0_dp/sqrt(9.0_dp*d)
    do
      z=rng_normal()
      if(z<=-1.0_dp/c)cycle
      x=d*(1.0_dp+c*z)**3
      call random_number(u)
      if(u<1.0_dp-0.0331_dp*z**4)exit
      if(log(u)<0.5_dp*z*z+d*(1.0_dp-x/d+log(x/d)))exit
    end do
    x=x*s
  end function rng_gamma

  real(dp) function rng_beta(a,b) result(x)
    real(dp),intent(in)::a,b
    real(dp)::g1,g2
    g1=rng_gamma(a);g2=rng_gamma(b);x=g1/(g1+g2)
  end function rng_beta

  real(dp) function rng_chisq(df) result(x)
    real(dp),intent(in)::df
    x=rng_gamma(0.5_dp*df,2.0_dp)
  end function rng_chisq

  integer function rng_poisson(lambda) result(k)
    real(dp),intent(in)::lambda
    real(dp)::u
    call random_number(u);k=poisson_quantile(u,lambda)
  end function rng_poisson

  integer function rng_binomial(n,p) result(k)
    integer,intent(in)::n
    real(dp),intent(in)::p
    real(dp)::u
    call random_number(u);k=binom_quantile(u,n,p)
  end function rng_binomial

  recursive real(dp) function adaptive_simpson(f,a,b,tol,max_depth) result(v)
    procedure(scalar_function)::f
    real(dp),intent(in)::a,b,tol
    integer,intent(in),optional::max_depth
    real(dp)::fa,fb,fc,s
    integer::depth
    depth=20;if(present(max_depth))depth=max_depth
    fa=f(a);fb=f(b);fc=f(0.5_dp*(a+b));s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    v=adaptive_simpson_rec(f,a,b,fa,fb,fc,s,tol,depth)
  end function adaptive_simpson

  recursive real(dp) function adaptive_simpson_rec(f,a,b,fa,fb,fc,s,tol,depth) result(v)
    procedure(scalar_function)::f
    real(dp),intent(in)::a,b,fa,fb,fc,s,tol
    integer,intent(in)::depth
    real(dp)::c,d,e,fd,fe,s1,s2
    c=0.5_dp*(a+b);d=0.5_dp*(a+c);e=0.5_dp*(c+b)
    fd=f(d);fe=f(e)
    s1=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
    s2=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
    if(depth<=0 .or. abs(s1+s2-s)<=15.0_dp*tol)then
      v=s1+s2+(s1+s2-s)/15.0_dp
    else
      v=adaptive_simpson_rec(f,a,c,fa,fc,fd,s1,0.5_dp*tol,depth-1)+ &
        adaptive_simpson_rec(f,c,b,fc,fb,fe,s2,0.5_dp*tol,depth-1)
    end if
  end function adaptive_simpson_rec

  real(dp) function bisect_root(f,a,b,tol,max_iter) result(x)
    procedure(scalar_function)::f
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::max_iter
    real(dp)::lo,hi,mid,fl,fm,t
    integer::i,n
    lo=a;hi=b;t=1.0e-10_dp;if(present(tol))t=tol
    n=200;if(present(max_iter))n=max_iter
    fl=f(lo)
    do i=1,n
      mid=0.5_dp*(lo+hi);fm=f(mid)
      if(fl*fm<=0.0_dp)then;hi=mid;else;lo=mid;fl=fm;end if
      if(abs(hi-lo)<=t*max(1.0_dp,abs(mid)))exit
    end do
    x=0.5_dp*(lo+hi)
  end function bisect_root

  subroutine solve_linear(a,b,x,info)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),intent(out)::x(:)
    integer,intent(out),optional::info
    real(dp),allocatable::aa(:,:),bb(:)
    real(dp)::factor,tmp
    integer::n,i,j,k,piv,istat
    n=size(b);allocate(aa(n,n),bb(n));aa=a;bb=b;istat=0
    do k=1,n-1
      piv=k
      do i=k+1,n
        if(abs(aa(i,k))>abs(aa(piv,k)))piv=i
      end do
      if(abs(aa(piv,k))<1.0e-14_dp)then;istat=1;exit;end if
      if(piv/=k)then
        do j=k,n;tmp=aa(k,j);aa(k,j)=aa(piv,j);aa(piv,j)=tmp;end do
        tmp=bb(k);bb(k)=bb(piv);bb(piv)=tmp
      end if
      do i=k+1,n
        factor=aa(i,k)/aa(k,k);aa(i,k:n)=aa(i,k:n)-factor*aa(k,k:n);bb(i)=bb(i)-factor*bb(k)
      end do
    end do
    if(istat==0 .and. abs(aa(n,n))<1.0e-14_dp)istat=1
    if(istat==0)then
      x(n)=bb(n)/aa(n,n)
      do i=n-1,1,-1;x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i);end do
    else
      x=0.0_dp
    end if
    if(present(info))info=istat
  end subroutine solve_linear

  subroutine invert_matrix(a,ainv,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::ainv(:,:)
    integer,intent(out),optional::info
    real(dp),allocatable::e(:),x(:)
    integer::n,j,istat,istatj
    n=size(a,1);allocate(e(n),x(n));ainv=0.0_dp;istat=0
    do j=1,n
      e=0.0_dp;e(j)=1.0_dp;call solve_linear(a,e,x,istatj)
      if(istatj/=0)then;istat=istatj;exit;end if
      ainv(:,j)=x
    end do
    if(present(info))info=istat
  end subroutine invert_matrix

  subroutine cholesky_lower(a,l,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(:,:)
    integer,intent(out)::info
    integer::n,i,j,k
    real(dp)::s
    n=size(a,1);l=0.0_dp;info=0
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1;s=s-l(i,k)*l(j,k);end do
        if(i==j)then
          if(s<=0.0_dp)then;info=i;return;end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine solve_spd(a,b,x,info)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),intent(out)::x(:)
    integer,intent(out)::info
    real(dp),allocatable::l(:,:),y(:)
    integer::n,i
    n=size(b);allocate(l(n,n),y(n));call cholesky_lower(a,l,info);if(info/=0)then;x=0.0_dp;return;end if
    do i=1,n;y(i)=(b(i)-dot_product(l(i,1:i-1),y(1:i-1)))/l(i,i);end do
    do i=n,1,-1;x(i)=(y(i)-dot_product(l(i+1:n,i),x(i+1:n)))/l(i,i);end do
  end subroutine solve_spd

  subroutine invert_spd(a,ainv,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::ainv(:,:)
    integer,intent(out)::info
    real(dp),allocatable::e(:),x(:)
    integer::n,j,istat
    n=size(a,1);allocate(e(n),x(n));ainv=0.0_dp;info=0
    do j=1,n
      e=0.0_dp;e(j)=1.0_dp;call solve_spd(a,e,x,istat)
      if(istat/=0)then;info=istat;return;end if
      ainv(:,j)=x
    end do
  end subroutine invert_spd

  subroutine mvn_random(mean,cov,x,info)
    real(dp),intent(in)::mean(:),cov(:,:)
    real(dp),intent(out)::x(:)
    integer,intent(out),optional::info
    real(dp),allocatable::l(:,:),z(:)
    integer::n,i,istat
    n=size(mean);allocate(l(n,n),z(n));call cholesky_lower(cov,l,istat)
    if(istat/=0)then;x=mean;if(present(info))info=istat;return;end if
    do i=1,n;z(i)=rng_normal();end do
    x=mean+matmul(l,z);if(present(info))info=0
  end subroutine mvn_random

end module tolerance_math
