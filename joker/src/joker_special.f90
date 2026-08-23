
module joker_special
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  integer, parameter, public :: dp = real64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  public :: normal_cdf, normal_quantile, reg_gamma_p, reg_beta, beta_quantile
  public :: digamma_j, trigamma_j, idigamma, rng_normal, rng_gamma, rng_beta
  public :: log_beta, median_real, variance_real
contains
  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: lo, hi, mid
    integer :: it
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp); return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp); return
    end if
    lo=-10.0_dp; hi=10.0_dp
    do it=1,100
      mid=0.5_dp*(lo+hi)
      if (normal_cdf(mid) < p) then
        lo=mid
      else
        hi=mid
      end if
    end do
    x=0.5_dp*(lo+hi)
  end function

  pure real(dp) function log_beta(a,b) result(v)
    real(dp), intent(in) :: a,b
    v = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
  end function

  pure real(dp) function reg_gamma_p(a,x) result(p)
    real(dp), intent(in) :: a,x
    integer, parameter :: itmax=500
    real(dp), parameter :: eps=1.0e-14_dp, fpmin=1.0e-300_dp
    real(dp) :: sumv, del, ap, b, c, d, h, an
    integer :: n
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan); return
    end if
    if (x <= 0.0_dp) then
      p=0.0_dp; return
    end if
    if (x < a+1.0_dp) then
      ap=a; sumv=1.0_dp/a; del=sumv
      do n=1,itmax
        ap=ap+1.0_dp; del=del*x/ap; sumv=sumv+del
        if (abs(del) < abs(sumv)*eps) exit
      end do
      p=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/b; h=d
      do n=1,itmax
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b; if (abs(d)<fpmin) d=fpmin
        c=b+an/c; if (abs(c)<fpmin) c=fpmin
        d=1.0_dp/d; del=d*c; h=h*del
        if (abs(del-1.0_dp)<eps) exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function

  pure real(dp) function beta_cf(a,b,x) result(h)
    real(dp), intent(in)::a,b,x
    integer, parameter :: maxit=500
    real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
    integer :: m,m2
    real(dp)::aa,c,d,del,qab,qam,qap
    qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
    c=1.0_dp; d=1.0_dp-qab*x/qap; if(abs(d)<fpmin)d=fpmin
    d=1.0_dp/d; h=d
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
      if(abs(del-1.0_dp)<eps)exit
    end do
  end function

  pure real(dp) function reg_beta(x,a,b) result(p)
    real(dp), intent(in)::x,a,b
    real(dp)::bt
    if(a<=0 .or. b<=0) then
      p=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    if(x<=0) then; p=0; return; end if
    if(x>=1) then; p=1; return; end if
    bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1-x))
    if(x < (a+1)/(a+b+2)) then
      p=bt*beta_cf(a,b,x)/a
    else
      p=1-bt*beta_cf(b,a,1-x)/b
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function

  pure real(dp) function beta_quantile(prob,a,b) result(x)
    real(dp), intent(in)::prob,a,b
    real(dp)::lo,hi,mid
    integer::it
    if(prob<=0)then;x=0;return;end if
    if(prob>=1)then;x=1;return;end if
    lo=0;hi=1
    do it=1,120
      mid=0.5_dp*(lo+hi)
      if(reg_beta(mid,a,b)<prob)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function

  pure real(dp) function digamma_j(x0) result(y)
    real(dp), intent(in)::x0
    real(dp)::x,r
    x=x0; y=0
    if(x<=0) then;y=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    do while(x<8); y=y-1/x; x=x+1; end do
    r=1/x
    y=y+log(x)-0.5_dp*r-r*r*(1.0_dp/12-r*r*(1.0_dp/120-r*r/252))
  end function

  pure real(dp) function trigamma_j(x0) result(y)
    real(dp), intent(in)::x0
    real(dp)::x,r
    x=x0;y=0
    if(x<=0) then;y=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    do while(x<8); y=y+1/(x*x);x=x+1;end do
    r=1/x
    y=y+r+0.5_dp*r*r+r**3/6-r**5/30+r**7/42
  end function

  pure real(dp) function idigamma(v) result(x)
    real(dp), intent(in)::v
    integer::it
    if(v>=-2.22_dp)then
      x=exp(v)+0.5_dp
    else
      x=-1.0_dp/(v+0.5772156649015329_dp)
    end if
    do it=1,20
      x=x-(digamma_j(x)-v)/trigamma_j(x)
      if(x<=0)x=1.0e-8_dp
    end do
  end function

  real(dp) function rng_normal() result(z)
    real(dp)::u1,u2
    call random_number(u1);call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2*log(u1))*cos(2*pi*u2)
  end function

  recursive real(dp) function rng_gamma(shape,scale) result(x)
    real(dp), intent(in)::shape,scale
    real(dp)::d,c,z,u
    if(shape<=0 .or. scale<=0) then
      x=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    if(shape<1) then
      call random_number(u)
      x=rng_gamma(shape+1,scale)*u**(1/shape)
      return
    end if
    d=shape-1.0_dp/3; c=1/sqrt(9*d)
    do
      z=rng_normal()
      if(1+c*z<=0)cycle
      x=(1+c*z)**3
      call random_number(u)
      if(u < 1-0.0331_dp*z**4)exit
      if(log(u) < 0.5_dp*z*z+d*(1-x+log(x)))exit
    end do
    x=scale*d*x
  end function

  real(dp) function rng_beta(a,b) result(x)
    real(dp),intent(in)::a,b
    real(dp)::u,v
    u=rng_gamma(a,1.0_dp);v=rng_gamma(b,1.0_dp);x=u/(u+v)
  end function

  pure real(dp) function median_real(x) result(m)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    real(dp)::tmp
    integer::i,j,n
    n=size(x);allocate(y(n));y=x
    do i=2,n
      tmp=y(i);j=i-1
      do while(j>=1)
        if(y(j)<=tmp)exit
        y(j+1)=y(j);j=j-1
      end do
      y(j+1)=tmp
    end do
    if(mod(n,2)==1) then;m=y((n+1)/2);else;m=0.5_dp*(y(n/2)+y(n/2+1));end if
  end function

  pure real(dp) function variance_real(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::m
    if(size(x)<2)then;v=0;return;end if
    m=sum(x)/size(x);v=sum((x-m)**2)/(size(x)-1)
  end function
end module joker_special
