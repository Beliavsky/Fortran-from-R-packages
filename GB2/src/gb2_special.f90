! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_special
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use gb2_kinds, only : dp, pi
  implicit none
  private
  public :: log_beta, beta_fn, digamma_fn, trigamma_fn, polygamma2_fn, polygamma3_fn
  public :: reg_incomplete_beta, beta_quantile, gamma_cdf, gamma_quantile
  public :: random_gamma, random_beta, adaptive_integral, hypergeo3f2_1, quiet_nan

  abstract interface
    function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface
contains
  pure real(dp) function quiet_nan() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan

  pure real(dp) function log_beta(a,b) result(v)
    real(dp), intent(in) :: a,b
    if(a<=0.0_dp .or. b<=0.0_dp) then
      v=quiet_nan()
    else
      v=log_gamma(a)+log_gamma(b)-log_gamma(a+b)
    end if
  end function log_beta

  pure real(dp) function beta_fn(a,b) result(v)
    real(dp), intent(in) :: a,b
    v=exp(log_beta(a,b))
  end function beta_fn

  pure real(dp) function digamma_fn(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: y, r, r2
    if(x<=0.0_dp) then
      v=quiet_nan()
      return
    end if
    y=x
    v=0.0_dp
    do while(y<10.0_dp)
      v=v-1.0_dp/y
      y=y+1.0_dp
    end do
    r=1.0_dp/y
    r2=r*r
    v=v+log(y)-0.5_dp*r-r2*(1.0_dp/12.0_dp-r2*(1.0_dp/120.0_dp-r2*(1.0_dp/252.0_dp-r2*(1.0_dp/240.0_dp-r2*5.0_dp/660.0_dp))))
  end function digamma_fn

  pure real(dp) function trigamma_fn(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: y,r,r2
    if(x<=0.0_dp) then
      v=quiet_nan()
      return
    end if
    y=x
    v=0.0_dp
    do while(y<10.0_dp)
      v=v+1.0_dp/(y*y)
      y=y+1.0_dp
    end do
    r=1.0_dp/y
    r2=r*r
    v=v+r+0.5_dp*r2+r*r2/6.0_dp-r*r2*r2/30.0_dp+r*r2*r2*r2/42.0_dp-r*r2*r2*r2*r2/30.0_dp+5.0_dp*r*r2**5/66.0_dp
  end function trigamma_fn

  pure real(dp) function polygamma2_fn(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: y,r
    if(x<=0.0_dp) then
      v=quiet_nan()
      return
    end if
    y=x
    v=0.0_dp
    do while(y<10.0_dp)
      v=v-2.0_dp/y**3
      y=y+1.0_dp
    end do
    r=1.0_dp/y
    v=v-r**2-r**3-0.5_dp*r**4+r**6/6.0_dp-r**8/6.0_dp+3.0_dp*r**10/10.0_dp-5.0_dp*r**12/6.0_dp
  end function polygamma2_fn

  pure real(dp) function polygamma3_fn(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: y,r
    if(x<=0.0_dp) then
      v=quiet_nan()
      return
    end if
    y=x
    v=0.0_dp
    do while(y<10.0_dp)
      v=v+6.0_dp/y**4
      y=y+1.0_dp
    end do
    r=1.0_dp/y
    v=v+2.0_dp*r**3+3.0_dp*r**4+2.0_dp*r**5-r**7+4.0_dp*r**9/3.0_dp-3.0_dp*r**11+10.0_dp*r**13
  end function polygamma3_fn

  pure real(dp) function betacf(a,b,x) result(h)
    real(dp), intent(in) :: a,b,x
    integer, parameter :: maxit=300
    real(dp), parameter :: fpmin=1.0e-300_dp, eps=3.0e-15_dp
    real(dp) :: qab,qap,qam,c,d,del,aa
    integer :: m,m2
    qab=a+b
    qap=a+1.0_dp
    qam=a-1.0_dp
    c=1.0_dp
    d=1.0_dp-qab*x/qap
    if(abs(d)<fpmin) d=fpmin
    d=1.0_dp/d
    h=d
    do m=1,maxit
      m2=2*m
      aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d=1.0_dp+aa*d
      if(abs(d)<fpmin) d=fpmin
      c=1.0_dp+aa/c
      if(abs(c)<fpmin) c=fpmin
      d=1.0_dp/d
      h=h*d*c
      aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
      d=1.0_dp+aa*d
      if(abs(d)<fpmin) d=fpmin
      c=1.0_dp+aa/c
      if(abs(c)<fpmin) c=fpmin
      d=1.0_dp/d
      del=d*c
      h=h*del
      if(abs(del-1.0_dp)<=eps) exit
    end do
  end function betacf

  pure real(dp) function reg_incomplete_beta(x,a,b) result(v)
    real(dp), intent(in) :: x,a,b
    real(dp) :: bt, xx
    if(a<=0.0_dp .or. b<=0.0_dp) then
      v=quiet_nan()
      return
    end if
    if(x<=0.0_dp) then
    v=0.0_dp
    return
    end if
    if(x>=1.0_dp) then
    v=1.0_dp
    return
    end if
    xx=min(1.0_dp,max(0.0_dp,x))
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(xx)+b*log(1.0_dp-xx))
    if(xx<(a+1.0_dp)/(a+b+2.0_dp)) then
      v=bt*betacf(a,b,xx)/a
    else
      v=1.0_dp-bt*betacf(b,a,1.0_dp-xx)/b
    end if
    v=min(1.0_dp,max(0.0_dp,v))
  end function reg_incomplete_beta

  real(dp) function beta_quantile(p,a,b) result(x)
    real(dp), intent(in) :: p,a,b
    real(dp) :: lo,hi,f,d,trial
    integer :: it
    if(a<=0.0_dp .or. b<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then
      x=quiet_nan()
      return
    end if
    if(p<=0.0_dp) then
    x=0.0_dp
    return
    end if
    if(p>=1.0_dp) then
    x=1.0_dp
    return
    end if
    lo=0.0_dp
    hi=1.0_dp
    x=a/(a+b)
    do it=1,100
      f=reg_incomplete_beta(x,a,b)-p
      if(f>0.0_dp) then
      hi=x
      else
      lo=x
      end if
      if(abs(f)<2.0e-14_dp .or. hi-lo<4.0e-15_dp) exit
      d=exp((a-1.0_dp)*log(max(x,tiny(1.0_dp)))+(b-1.0_dp)*log(max(1.0_dp-x,tiny(1.0_dp)))-log_beta(a,b))
      if(d>0.0_dp .and. ieee_is_finite(d)) then
        trial=x-f/d
      else
        trial=0.5_dp*(lo+hi)
      end if
      if(trial<=lo .or. trial>=hi .or. .not.ieee_is_finite(trial)) trial=0.5_dp*(lo+hi)
      x=trial
    end do
    x=min(1.0_dp,max(0.0_dp,x))
  end function beta_quantile

  pure real(dp) function gamma_cdf(x,a) result(p)
    real(dp), intent(in) :: x,a
    integer, parameter :: itmax=500
    real(dp), parameter :: eps=3.0e-15_dp, fpmin=1.0e-300_dp
    real(dp) :: ap,sumv,del,b,c,d,h,an,q
    integer :: n
    if(a<=0.0_dp) then
    p=quiet_nan()
    return
    end if
    if(x<=0.0_dp) then
    p=0.0_dp
    return
    end if
    if(x<a+1.0_dp) then
      ap=a
      sumv=1.0_dp/a
      del=sumv
      do n=1,itmax
        ap=ap+1.0_dp
        del=del*x/ap
        sumv=sumv+del
        if(abs(del)<abs(sumv)*eps) exit
      end do
      p=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a
      c=1.0_dp/fpmin
      d=1.0_dp/b
      h=d
      do n=1,itmax
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b
        if(abs(d)<fpmin) d=fpmin
        c=b+an/c
        if(abs(c)<fpmin) c=fpmin
        d=1.0_dp/d
        del=d*c
        h=h*del
        if(abs(del-1.0_dp)<eps) exit
      end do
      q=exp(-x+a*log(x)-log_gamma(a))*h
      p=1.0_dp-q
    end if
    p=min(1.0_dp,max(0.0_dp,p))
  end function gamma_cdf

  real(dp) function gamma_quantile(p,a) result(x)
    real(dp), intent(in) :: p,a
    real(dp) :: lo,hi,f,d,trial
    integer :: it
    if(a<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then
    x=quiet_nan()
    return
    end if
    if(p<=0.0_dp) then
    x=0.0_dp
    return
    end if
    if(p>=1.0_dp) then
    x=huge(1.0_dp)
    return
    end if
    lo=0.0_dp
    hi=max(1.0_dp,a)
    do while(gamma_cdf(hi,a)<p .and. hi<huge(1.0_dp)/4.0_dp)
      hi=2.0_dp*hi
    end do
    x=min(hi,max(lo,a))
    do it=1,120
      f=gamma_cdf(x,a)-p
      if(f>0.0_dp) then
      hi=x
      else
      lo=x
      end if
      if(abs(f)<2.0e-14_dp .or. hi-lo<1.0e-13_dp*max(1.0_dp,x)) exit
      if(x>0.0_dp) then
        d=exp((a-1.0_dp)*log(x)-x-log_gamma(a))
      else
        d=0.0_dp
      end if
      if(d>0.0_dp .and. ieee_is_finite(d)) then
      trial=x-f/d
      else
      trial=0.5_dp*(lo+hi)
      end if
      if(trial<=lo .or. trial>=hi .or. .not.ieee_is_finite(trial)) trial=0.5_dp*(lo+hi)
      x=trial
    end do
  end function gamma_quantile

  recursive subroutine random_gamma(shape,x)
    real(dp), intent(in) :: shape
    real(dp), intent(out) :: x
    real(dp) :: d,c,z,u,v,g
    if(shape<=0.0_dp) then
    x=quiet_nan()
    return
    end if
    if(shape<1.0_dp) then
      call random_gamma(shape+1.0_dp,g)
      call random_number(u)
      x=g*u**(1.0_dp/shape)
      return
    end if
    d=shape-1.0_dp/3.0_dp
    c=1.0_dp/sqrt(9.0_dp*d)
    do
      call normal_random(z)
      v=(1.0_dp+c*z)**3
      if(v<=0.0_dp) cycle
      call random_number(u)
      if(u<1.0_dp-0.0331_dp*z**4) exit
      if(log(max(u,tiny(1.0_dp)))<0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x=d*v
  end subroutine random_gamma

  subroutine normal_random(z)
    real(dp), intent(out) :: z
    real(dp) :: u1,u2
    call random_number(u1)
    call random_number(u2)
    z=sqrt(-2.0_dp*log(max(u1,tiny(1.0_dp))))*cos(2.0_dp*pi*u2)
  end subroutine normal_random

  subroutine random_beta(a,b,x)
    real(dp), intent(in) :: a,b
    real(dp), intent(out) :: x
    real(dp) :: g1,g2
    call random_gamma(a,g1)
    call random_gamma(b,g2)
    x=g1/(g1+g2)
  end subroutine random_beta

  function adaptive_integral(f,a,b,tol,max_depth) result(value)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_depth
    real(dp) :: value,err,target
    integer :: depth
    target=1.0e-10_dp
    if(present(tol)) target=max(tol,10.0_dp*epsilon(1.0_dp))
    depth=25
    if(present(max_depth)) depth=max(1,max_depth)
    if(abs(b-a)<=tiny(1.0_dp)) then
    value=0.0_dp
    return
    end if
    call adapt_gk(f,a,b,target,depth,value,err)
  end function adaptive_integral

  recursive subroutine adapt_gk(f,a,b,tol,depth,val,err)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b,tol
    integer, intent(in) :: depth
    real(dp), intent(out) :: val,err
    real(dp) :: v0,e0,v1,e1,v2,e2,m
    call gk15(f,a,b,v0,e0)
    if(depth<=0 .or. e0<=tol*max(1.0_dp,abs(v0))) then
      val=v0
      err=e0
      return
    end if
    m=0.5_dp*(a+b)
    call adapt_gk(f,a,m,0.5_dp*tol,depth-1,v1,e1)
    call adapt_gk(f,m,b,0.5_dp*tol,depth-1,v2,e2)
    val=v1+v2
    err=e1+e2
  end subroutine adapt_gk

  subroutine gk15(f,a,b,val,err)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a,b
    real(dp), intent(out) :: val,err
    real(dp), parameter :: xgk(8)=[0.9914553711208126_dp,0.9491079123427585_dp,0.8648644233597691_dp,0.7415311855993945_dp, &
      0.5860872354676911_dp,0.4058451513773972_dp,0.2077849550078985_dp,0.0_dp]
    real(dp), parameter :: wgk(8)=[0.02293532201052922_dp,0.06309209262997855_dp,0.1047900103222502_dp,0.1406532597155259_dp, &
      0.1690047266392679_dp,0.1903505780647854_dp,0.2044329400752989_dp,0.2094821410847278_dp]
    real(dp), parameter :: wg(4)=[0.1294849661688697_dp,0.2797053914892767_dp,0.3818300505051189_dp,0.4179591836734694_dp]
    real(dp) :: c,h,fc,resg,resk,f1,f2,absc
    integer :: j
    c=0.5_dp*(a+b)
    h=0.5_dp*(b-a)
    fc=f(c)
    resg=wg(4)*fc
    resk=wgk(8)*fc
    do j=1,7
      absc=h*xgk(j)
      f1=f(c-absc)
      f2=f(c+absc)
      resk=resk+wgk(j)*(f1+f2)
      select case(j)
      case(2)
      resg=resg+wg(1)*(f1+f2)
      case(4)
      resg=resg+wg(2)*(f1+f2)
      case(6)
      resg=resg+wg(3)*(f1+f2)
      end select
    end do
    val=resk*h
    err=abs((resk-resg)*h)
  end subroutine gk15

  real(dp) function hypergeo3f2_1(u,l,tol,maxiter,converged) result(v)
    real(dp), intent(in) :: u(3),l(2)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    logical, intent(out), optional :: converged
    real(dp) :: term,t
    integer :: n,nmax
    logical :: ok
    t=1.0e-12_dp
    if(present(tol)) t=max(tol,10.0_dp*epsilon(1.0_dp))
    nmax=100000
    if(present(maxiter)) nmax=max(1,maxiter)
    term=1.0_dp
    v=1.0_dp
    ok=.false.
    do n=1,nmax
      term=term*((u(1)+n-1.0_dp)*(u(2)+n-1.0_dp)*(u(3)+n-1.0_dp))/ &
        ((l(1)+n-1.0_dp)*(l(2)+n-1.0_dp)*real(n,dp))
      v=v+term
      if(.not.ieee_is_finite(v)) exit
      if(abs(term)<=t*max(1.0_dp,abs(v))) then
      ok=.true.
      exit
      end if
    end do
    if(present(converged)) converged=ok
  end function hypergeo3f2_1
end module gb2_special
