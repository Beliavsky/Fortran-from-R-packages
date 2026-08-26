! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_math
  use waveslim_kinds, only : dp, pi, i8
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use r_descriptive, only : r_mean, r_variance
  use r_distributions, only : r_pchisq, r_pnorm, r_qnorm
  use r_quantiles, only : r_median, r_quantile_type7
  use r_robust, only : r_mad
  implicit none
  private
  public :: mean_value, variance_value, median_value, mad_value, quantile_type7
  public :: normal_cdf, normal_quantile, chi_square_cdf, fft_complex, next_power_two
  public :: seed_rng, random_normal, circular_shift, linear_interp, trapz
  public :: autocovariance, autocorrelation, cross_correlation, lower_string
contains
  pure function mean_value(x) result(ans)
    real(dp), intent(in) :: x(:)
    real(dp) :: ans
    real(dp), allocatable :: finite_x(:)
    finite_x=pack(x,ieee_is_finite(x))
    if(size(finite_x)==0) then
      ans=0.0_dp
    else
      ans=r_mean(finite_x)
    end if
  end function mean_value

  pure function variance_value(x, sample) result(ans)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: sample
    real(dp) :: ans
    real(dp), allocatable :: finite_x(:)
    integer :: den
    logical :: smp
    smp=.true.
    if(present(sample)) smp=sample
    finite_x=pack(x,ieee_is_finite(x))
    den=merge(1,0,smp)
    if(size(finite_x)>den) then
      ans=r_variance(finite_x,ddof=den)
    else
      ans=0.0_dp
    end if
  end function variance_value

  function median_value(x) result(ans)
    real(dp), intent(in) :: x(:)
    real(dp) :: ans
    real(dp), allocatable :: finite_x(:)
    finite_x=pack(x,ieee_is_finite(x))
    if(size(finite_x)==0) then
      ans=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      ans=r_median(finite_x)
    end if
  end function median_value

  function mad_value(x, center, constant) result(ans)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: center, constant
    real(dp) :: ans,c,k
    real(dp), allocatable :: finite_x(:)
    finite_x=pack(x,ieee_is_finite(x))
    c=r_median(finite_x)
    if(present(center)) c=center
    k=1.482602218505602_dp
    if(present(constant)) k=constant
    if(.not.ieee_is_finite(c)) then
      ans=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      ans=r_mad(finite_x,center=c,constant=k)
    end if
  end function mad_value

  function quantile_type7(x,p) result(ans)
    real(dp), intent(in) :: x(:),p
    real(dp) :: ans
    real(dp),allocatable::finite_x(:)
    finite_x=pack(x,ieee_is_finite(x))
    if(size(finite_x)==0 .or. p<0.0_dp .or. p>1.0_dp) then
      ans=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      ans=r_quantile_type7(finite_x,p)
    end if
  end function quantile_type7

  pure elemental function normal_cdf(x) result(p)
    real(dp),intent(in)::x
    real(dp)::p
    p=r_pnorm(x)
  end function normal_cdf

  pure elemental function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x
    if(p<=0.0_dp) then
      x=-huge(1.0_dp)
    else if(p>=1.0_dp) then
      x=huge(1.0_dp)
    else
      x=r_qnorm(p)
    end if
  end function normal_quantile

  function chi_square_cdf(x,df) result(p)
    real(dp),intent(in)::x,df
    real(dp)::p
    if(x<=0.0_dp) then
    p=0.0_dp
    else
    p=r_pchisq(x,df)
    end if
  end function chi_square_cdf

  subroutine seed_rng(seed)
    integer(i8),intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
    put(i)=int(mod(seed+104729_i8*int(i,i8),int(huge(1),i8)-1_i8))+1
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  subroutine random_normal(x)
    real(dp),intent(out)::x(:)
    real(dp)::u1,u2
    integer::i
    i=1
    do while(i<=size(x))
      call random_number(u1)
      call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      x(i)=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
      if(i+1<=size(x)) x(i+1)=sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
      i=i+2
    end do
  end subroutine random_normal

  function circular_shift(x,m) result(y)
    real(dp),intent(in)::x(:)
    integer,intent(in)::m
    real(dp),allocatable::y(:)
    integer::i,n,k
    n=size(x)
    allocate(y(n))
    do i=1,n
    k=modulo(i-1+m,n)+1
    y(i)=x(k)
    end do
  end function circular_shift

  pure function linear_interp(x,y,x0) result(y0)
    real(dp),intent(in)::x(:),y(:),x0
    real(dp)::y0,t
    integer::i,n
    n=size(x)
    if(x0<=x(1)) then
    y0=y(1)
    return
    else if(x0>=x(n)) then
    y0=y(n)
    return
    end if
    do i=1,n-1
      if(x0<=x(i+1)) then
      t=(x0-x(i))/(x(i+1)-x(i))
      y0=(1-t)*y(i)+t*y(i+1)
      return
      end if
    end do
    y0=y(n)
  end function linear_interp

  pure function trapz(x,y) result(ans)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::ans
    integer::i
    ans=0.0_dp
    do i=1,size(x)-1
    ans=ans+0.5_dp*(x(i+1)-x(i))*(y(i+1)+y(i))
    end do
  end function trapz

  function autocovariance(x,lag,max_unbiased) result(v)
    real(dp),intent(in)::x(:)
    integer,intent(in)::lag
    logical,intent(in),optional::max_unbiased
    real(dp)::v,m
    integer::i,n,den
    logical::ub
    ub=.true.
    if(present(max_unbiased))ub=max_unbiased
    n=size(x)
    m=mean_value(x)
    v=0.0_dp
    do i=1,n-lag
    v=v+(x(i)-m)*(x(i+lag)-m)
    end do
    den=n-lag
    if(.not.ub)den=n
    if(den>0)v=v/real(den,dp)
  end function autocovariance

  function autocorrelation(x,lag) result(r)
    real(dp),intent(in)::x(:)
    integer,intent(in)::lag
    real(dp)::r,v0
    v0=autocovariance(x,0)
    if(v0>0.0_dp) then
    r=autocovariance(x,lag)/v0
    else
    r=0.0_dp
    end if
  end function autocorrelation

  function cross_correlation(x,y,lag) result(r)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in)::lag
    real(dp)::r,mx,my,sx,sy,c
    integer::i,n,nuse
    n=min(size(x),size(y))
    mx=mean_value(x(1:n))
    my=mean_value(y(1:n))
    c=0.0_dp
    nuse=0
    if(lag>=0)then
      do i=1,n-lag
      c=c+(x(i)-mx)*(y(i+lag)-my)
      nuse=nuse+1
      end do
    else
      do i=1,n+lag
      c=c+(x(i-lag)-mx)*(y(i)-my)
      nuse=nuse+1
      end do
    end if
    sx=sqrt(max(variance_value(x(1:n),.false.),0.0_dp))
    sy=sqrt(max(variance_value(y(1:n),.false.),0.0_dp))
    if(nuse>0 .and. sx*sy>0.0_dp)then
    r=c/real(nuse,dp)/(sx*sy)
    else
    r=0.0_dp
    end if
  end function cross_correlation

  subroutine fft_complex(a,inverse)
    complex(dp),intent(inout)::a(:)
    logical,intent(in),optional::inverse
    integer::n,i,j,m,mmax,istep
    complex(dp)::temp,w,wm
    logical::inv
    real(dp)::theta
    n=size(a)
    inv=.false.
    if(present(inverse))inv=inverse
    if(iand(n,n-1)/=0)error stop 'fft_complex requires power-of-two length'
    j=1
    do i=1,n
      if(j>i)then
      temp=a(j)
      a(j)=a(i)
      a(i)=temp
      end if
      m=n/2
      do while(m>=1 .and. j>m)
      j=j-m
      m=m/2
      end do
      j=j+m
    end do
    mmax=1
    do while(n>mmax)
      istep=2*mmax
      theta=merge(2.0_dp,-2.0_dp,inv)*pi/real(istep,dp)
      wm=cmplx(cos(theta),sin(theta),dp)
      w=(1.0_dp,0.0_dp)
      do m=1,mmax
        do i=m,n,istep
        j=i+mmax
        temp=w*a(j)
        a(j)=a(i)-temp
        a(i)=a(i)+temp
        end do
        w=w*wm
      end do
      mmax=istep
    end do
    if(inv)a=a/real(n,dp)
  end subroutine fft_complex

  integer pure function next_power_two(n) result(p)
    integer,intent(in)::n
    p=1
    do while(p<n)
    p=2*p
    end do
  end function next_power_two

  pure function lower_string(text) result(ans)
    character(len=*),intent(in)::text
    character(len=len(text))::ans
    integer::i,c
    do i=1,len(text)
    c=iachar(text(i:i))
    if(c>=65.and.c<=90)then
    ans(i:i)=achar(c+32)
    else
    ans(i:i)=text(i:i)
    end if
    end do
  end function lower_string
end module waveslim_math
