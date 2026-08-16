! SPDX-License-Identifier: GPL-2.0-only
module mbbefd_math
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use mbbefd_kinds, only : dp
  use fitdistrplus_math, only : regularized_beta, random_gamma, random_poisson, type7_quantile, &
    invert_matrix, numerical_hessian
  implicit none
  private

  public :: nan_dp, log_beta, beta_pdf, beta_cdf, beta_quantile, beta_random
  public :: ncbeta_pdf, ncbeta_cdf, ncbeta_quantile, ncbeta_random
  public :: probability_to_lower, gauss_legendre_integrate, digamma_dp
  public :: type7_quantile, invert_matrix, numerical_hessian

  abstract interface
    function integrand_callback(x, context) result(value)
      import dp
      real(dp), intent(in) :: x
      class(*), intent(in) :: context
      real(dp) :: value
    end function integrand_callback
  end interface

contains

  pure function nan_dp() result(x)
    real(dp) :: x
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure function log_beta(a,b) result(v)
    real(dp), intent(in) :: a,b
    real(dp) :: v
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      v = nan_dp()
    else
      v = log_gamma(a) + log_gamma(b) - log_gamma(a+b)
    end if
  end function log_beta

  pure function beta_pdf(x,a,b,log_value) result(v)
    real(dp), intent(in) :: x,a,b
    logical, intent(in), optional :: log_value
    real(dp) :: v,lv
    logical :: want_log
    want_log = .false.; if (present(log_value)) want_log = log_value
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      v = nan_dp(); return
    end if
    if (x < 0.0_dp .or. x > 1.0_dp) then
      v = merge(-huge(1.0_dp),0.0_dp,want_log); return
    end if
    if (x == 0.0_dp) then
      if (a < 1.0_dp) then
        v = huge(1.0_dp)
      else if (a == 1.0_dp) then
        lv = -log_beta(a,b); v = merge(lv,exp(lv),want_log)
      else
        v = merge(-huge(1.0_dp),0.0_dp,want_log)
      end if
      return
    end if
    if (x == 1.0_dp) then
      if (b < 1.0_dp) then
        v = huge(1.0_dp)
      else if (b == 1.0_dp) then
        lv = -log_beta(a,b); v = merge(lv,exp(lv),want_log)
      else
        v = merge(-huge(1.0_dp),0.0_dp,want_log)
      end if
      return
    end if
    lv = (a-1.0_dp)*log(x) + (b-1.0_dp)*log(1.0_dp-x) - log_beta(a,b)
    v = merge(lv,exp(lv),want_log)
  end function beta_pdf

  pure function beta_cdf(x,a,b) result(v)
    real(dp), intent(in) :: x,a,b
    real(dp) :: v
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      v = nan_dp()
    else if (x <= 0.0_dp) then
      v = 0.0_dp
    else if (x >= 1.0_dp) then
      v = 1.0_dp
    else
      v = regularized_beta(x,a,b)
    end if
  end function beta_cdf

  pure function beta_quantile(p,a,b) result(x)
    real(dp), intent(in) :: p,a,b
    real(dp) :: x,lo,hi,mid,f
    integer :: it
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
      x = nan_dp(); return
    end if
    if (p == 0.0_dp) then; x=0.0_dp; return; end if
    if (p == 1.0_dp) then; x=1.0_dp; return; end if
    lo=0.0_dp; hi=1.0_dp
    do it=1,120
      mid=0.5_dp*(lo+hi); f=regularized_beta(mid,a,b)
      if (f < p) then; lo=mid; else; hi=mid; end if
      if (hi-lo <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(mid))) exit
    end do
    x=0.5_dp*(lo+hi)
  end function beta_quantile

  function beta_random(a,b) result(x)
    real(dp), intent(in) :: a,b
    real(dp) :: x,g1,g2
    if (a <= 0.0_dp .or. b <= 0.0_dp) then; x=nan_dp(); return; end if
    g1=random_gamma(a,1.0_dp); g2=random_gamma(b,1.0_dp)
    if (g1+g2 > 0.0_dp) then; x=g1/(g1+g2); else; x=0.5_dp; end if
  end function beta_random

  function ncbeta_pdf(x,a,b,ncp,log_value) result(v)
    real(dp), intent(in) :: x,a,b,ncp
    logical, intent(in), optional :: log_value
    real(dp) :: v,w,sumv,lambda,term
    integer :: k
    logical :: want_log
    want_log=.false.; if(present(log_value)) want_log=log_value
    if (ncp < 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then; v=nan_dp(); return; end if
    if (ncp == 0.0_dp) then; v=beta_pdf(x,a,b,want_log); return; end if
    if (x < 0.0_dp .or. x > 1.0_dp) then; v=merge(-huge(1.0_dp),0.0_dp,want_log); return; end if
    lambda=0.5_dp*ncp
    if (lambda > 700.0_dp) then
      ! The upstream package only uses ncp=0 in its exposure/moment layer.
      ! Fall back to a direct simulation-style mode-centered mixture is not needed here.
      v=nan_dp(); return
    end if
    w=exp(-lambda); sumv=0.0_dp
    do k=0,100000
      term=beta_pdf(x,a+real(k,dp),b)
      sumv=sumv+w*term
      if (k > lambda .and. w*max(1.0_dp,abs(term)) < 1.0e-15_dp*max(1.0_dp,abs(sumv))) exit
      w=w*lambda/real(k+1,dp)
      if (w == 0.0_dp) exit
    end do
    if (want_log) then
      if(sumv>0.0_dp) then; v=log(sumv); else; v=-huge(1.0_dp); end if
    else
      v=sumv
    end if
  end function ncbeta_pdf

  function ncbeta_cdf(x,a,b,ncp) result(v)
    real(dp), intent(in) :: x,a,b,ncp
    real(dp) :: v,w,lambda,term
    integer :: k
    if (ncp < 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then; v=nan_dp(); return; end if
    if (x <= 0.0_dp) then; v=0.0_dp; return; end if
    if (x >= 1.0_dp) then; v=1.0_dp; return; end if
    if (ncp == 0.0_dp) then; v=beta_cdf(x,a,b); return; end if
    lambda=0.5_dp*ncp
    if (lambda > 700.0_dp) then; v=nan_dp(); return; end if
    w=exp(-lambda); v=0.0_dp
    do k=0,100000
      term=beta_cdf(x,a+real(k,dp),b)
      v=v+w*term
      if (k > lambda .and. w < 1.0e-15_dp) exit
      w=w*lambda/real(k+1,dp)
      if (w == 0.0_dp) exit
    end do
    v=min(1.0_dp,max(0.0_dp,v))
  end function ncbeta_cdf

  function ncbeta_quantile(p,a,b,ncp) result(x)
    real(dp), intent(in) :: p,a,b,ncp
    real(dp) :: x,lo,hi,mid
    integer :: it
    if (p < 0.0_dp .or. p > 1.0_dp .or. a<=0.0_dp .or. b<=0.0_dp .or. ncp<0.0_dp) then
      x=nan_dp(); return
    end if
    if(p==0.0_dp) then; x=0.0_dp; return; end if
    if(p==1.0_dp) then; x=1.0_dp; return; end if
    lo=0.0_dp;hi=1.0_dp
    do it=1,120
      mid=0.5_dp*(lo+hi)
      if(ncbeta_cdf(mid,a,b,ncp)<p) then; lo=mid; else; hi=mid; end if
    end do
    x=0.5_dp*(lo+hi)
  end function ncbeta_quantile

  function ncbeta_random(a,b,ncp) result(x)
    real(dp), intent(in) :: a,b,ncp
    real(dp) :: x,g1,g2,lambda
    integer :: k
    if(a<=0.0_dp .or. b<=0.0_dp .or. ncp<0.0_dp) then; x=nan_dp(); return; end if
    lambda=0.5_dp*ncp; k=random_poisson(lambda)
    g1=random_gamma(a+real(k,dp),1.0_dp); g2=random_gamma(b,1.0_dp)
    if(g1+g2>0.0_dp) then; x=g1/(g1+g2); else; x=0.5_dp; end if
  end function ncbeta_random

  pure function digamma_dp(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v,y,r,inv,inv2
    if (x <= 0.0_dp) then
      v=nan_dp(); return
    end if
    y=x; r=0.0_dp
    do while(y < 8.0_dp)
      r=r-1.0_dp/y; y=y+1.0_dp
    end do
    inv=1.0_dp/y; inv2=inv*inv
    v=r+log(y)-0.5_dp*inv-inv2*(1.0_dp/12.0_dp-inv2*(1.0_dp/120.0_dp- &
      inv2*(1.0_dp/252.0_dp-inv2*(1.0_dp/240.0_dp))))
  end function digamma_dp

  pure function probability_to_lower(p,lower_tail,log_p) result(v)
    real(dp), intent(in) :: p
    logical, intent(in), optional :: lower_tail,log_p
    real(dp) :: v
    logical :: low,lp
    low=.true.;lp=.false.;if(present(lower_tail))low=lower_tail;if(present(log_p))lp=log_p
    if(lp) then; v=exp(p); else; v=p; end if
    if(.not.low) v=1.0_dp-v
  end function probability_to_lower

  function gauss_legendre_integrate(fun,context,a,b,n) result(value)
    procedure(integrand_callback) :: fun
    class(*), intent(in) :: context
    real(dp), intent(in) :: a,b
    integer, intent(in), optional :: n
    real(dp) :: value,z,z1,p1,p2,p3,pp,xm,xl
    integer :: nn,m,i,j
    nn=96; if(present(n)) nn=max(8,n)
    m=(nn+1)/2; xm=0.5_dp*(b+a); xl=0.5_dp*(b-a); value=0.0_dp
    do i=1,m
      z=cos(acos(-1.0_dp)*(real(i,dp)-0.25_dp)/(real(nn,dp)+0.5_dp))
      do
        p1=1.0_dp;p2=0.0_dp
        do j=1,nn
          p3=p2;p2=p1;p1=((2.0_dp*real(j,dp)-1.0_dp)*z*p2-(real(j,dp)-1.0_dp)*p3)/real(j,dp)
        end do
        pp=real(nn,dp)*(z*p1-p2)/(z*z-1.0_dp)
        z1=z;z=z1-p1/pp
        if(abs(z-z1)<=4.0_dp*epsilon(1.0_dp))exit
      end do
      value=value+2.0_dp*xl/((1.0_dp-z*z)*pp*pp)*(fun(xm-xl*z,context)+fun(xm+xl*z,context))
    end do
  end function gauss_legendre_integrate

end module mbbefd_math
