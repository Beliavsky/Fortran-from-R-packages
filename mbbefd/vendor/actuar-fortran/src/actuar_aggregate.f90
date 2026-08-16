! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_aggregate
  use actuar_kinds, only : dp, i8
  use actuar_special, only : nan_dp, normal_cdf, normal_pdf, normal_quantile
  use actuar_rng, only : seed_rng
  use actuar_types, only : aggregate_distribution
  implicit none
  private

  abstract interface
    function cdf_callback(x) result(p)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: p
    end function cdf_callback
    function lev_callback(x) result(m)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: m
    end function lev_callback
    function integer_rng_callback() result(k)
      integer :: k
    end function integer_rng_callback
    function real_rng_callback() result(x)
      import dp
      real(dp) :: x
    end function real_rng_callback
  end interface

  public :: discretize_cdf, exact_compound, panjer_recursion
  public :: panjer_poisson, panjer_binomial, panjer_nbinomial
  public :: compound_simulation, aggregate_var, aggregate_cte
  public :: normal_aggregate_cdf, normal_power_cdf, convolve_probabilities

contains

  function discretize_cdf(cdf, from, to, step, method, lev) result(prob)
    procedure(cdf_callback) :: cdf
    real(dp), intent(in) :: from, to, step
    character(len=*), intent(in) :: method
    procedure(lev_callback), optional :: lev
    real(dp), allocatable :: prob(:)
    real(dp), allocatable :: x(:), f(:), e(:)
    integer :: n, i

    if (step <= 0.0_dp .or. to <= from) then
      allocate(prob(0)); return
    end if
    n = nint((to-from)/step)
    if (n < 1) then
      allocate(prob(0))
      return
    end if

    select case(trim(adjustl(method)))
    case('upper')
      allocate(prob(n),x(n+1),f(n+1))
      x = [(from+real(i,dp)*step,i=0,n)]
      do i=1,n+1; f(i)=cdf(x(i)); end do
      prob = f(2:)-f(:n)
    case('lower')
      allocate(prob(n+1),x(n+1),f(n+1))
      x = [(from+real(i,dp)*step,i=0,n)]
      do i=1,n+1; f(i)=cdf(x(i)); end do
      prob(1)=0.0_dp; prob(2:)=f(2:)-f(:n)
    case('rounding')
      allocate(prob(n),x(n+1),f(n+1))
      x(1)=from
      do i=2,n+1
        x(i)=from+(real(i-1,dp)-0.5_dp)*step
      end do
      do i=1,n+1; f(i)=cdf(x(i)); end do
      prob=f(2:)-f(:n)
    case('unbiased')
      if (.not. present(lev)) then
        allocate(prob(0)); return
      end if
      allocate(prob(n+1),x(n+1),e(n+1))
      prob = 0.0_dp
      x = 0.0_dp
      e = 0.0_dp
      x=[(from+real(i,dp)*step,i=0,n)]
      do i=1,n+1; e(i)=lev(x(i)); end do
      block
        real(dp) :: e_first, e_second
        e_first = e(1)
        e_second = e(2)
        prob(1) = -(e_second-e_first)/step + 1.0_dp - cdf(from)
      end block
      do i=2,n
        prob(i)=(2.0_dp*e(i)-e(i-1)-e(i+1))/step
      end do
      prob(n+1)=(e(n+1)-e(n))/step-1.0_dp+cdf(to)
    case default
      allocate(prob(0)); return
    end select
    where (abs(prob) < 1.0e-15_dp) prob=0.0_dp
  end function discretize_cdf

  function convolve_probabilities(a,b) result(c)
    real(dp), intent(in) :: a(:),b(:)
    real(dp), allocatable :: c(:)
    integer :: i,j
    allocate(c(size(a)+size(b)-1)); c=0.0_dp
    do i=1,size(a)
      do j=1,size(b)
        c(i+j-1)=c(i+j-1)+a(i)*b(j)
      end do
    end do
  end function convolve_probabilities

  function exact_compound(severity, frequency, scale) result(out)
    real(dp), intent(in) :: severity(:),frequency(:)
    real(dp), intent(in), optional :: scale
    type(aggregate_distribution) :: out
    real(dp), allocatable :: current(:),temp(:)
    real(dp) :: h
    integer :: nmax,total_len,n,i
    h=1.0_dp; if(present(scale)) h=scale
    nmax=size(frequency)-1
    total_len=nmax*(size(severity)-1)+1
    allocate(out%probability(total_len)); out%probability=0.0_dp
    out%probability(1)=frequency(1)
    allocate(current(1)); current=1.0_dp
    do n=1,nmax
      temp=convolve_probabilities(current,severity)
      call move_alloc(temp,current)
      out%probability(:size(current))=out%probability(:size(current))+frequency(n+1)*current
    end do
    allocate(out%support(total_len),out%cdf(total_len))
    out%support=[(real(i-1,dp)*h,i=1,total_len)]
    out%cdf(1)=out%probability(1)
    do i=2,total_len; out%cdf(i)=out%cdf(i-1)+out%probability(i); end do
    out%cdf=min(1.0_dp,out%cdf)
    out%ok=.true.; out%message='exact convolution'
  end function exact_compound

  function panjer_recursion(severity,a,b,fs0,p0,p1,tolerance,max_terms,scale) result(out)
    real(dp), intent(in) :: severity(:),a,b,fs0
    real(dp), intent(in), optional :: p0,p1,tolerance,scale
    integer, intent(in), optional :: max_terms
    type(aggregate_distribution) :: out
    real(dp), allocatable :: fs(:)
    real(dp) :: norm,cumul,tol,h,term,extra,pzero,pone
    integer :: maxn,x,k,m,i
    logical :: ab1
    maxn=500; if(present(max_terms)) maxn=max_terms
    tol=1.0_dp-1.0e-8_dp; if(present(tolerance)) tol=tolerance
    h=1.0_dp; if(present(scale)) h=scale
    allocate(fs(maxn+1)); fs=0.0_dp; fs(1)=fs0
    cumul=fs0; norm=1.0_dp-a*severity(1)
    if(abs(norm)<1.0e-14_dp) then
      out%ok=.false.; out%message='singular Panjer normalization'; return
    end if
    ab1=present(p0)
    pzero=0.0_dp; pone=0.0_dp
    if(ab1) pzero=p0
    if(present(p1)) pone=p1
    term=pone-(a+b)*pzero
    x=1
    do while(x<=maxn .and. cumul<tol)
      m=min(x,size(severity)-1)
      do k=1,m
        fs(x+1)=fs(x+1)+(a+b*real(k,dp)/real(x,dp))*severity(k+1)*fs(x-k+1)
      end do
      extra=0.0_dp
      if(ab1 .and. x<=size(severity)-1) extra=severity(x+1)*term
      fs(x+1)=(fs(x+1)+extra)/norm
      if(fs(x+1)<0.0_dp .and. abs(fs(x+1))<1.0e-13_dp) fs(x+1)=0.0_dp
      cumul=cumul+fs(x+1)
      x=x+1
    end do
    allocate(out%probability(x),out%support(x),out%cdf(x))
    out%probability=fs(:x)
    out%support=[(real(i-1,dp)*h,i=1,x)]
    out%cdf(1)=out%probability(1)
    do i=2,x; out%cdf(i)=out%cdf(i-1)+out%probability(i); end do
    out%cdf=min(1.0_dp,out%cdf)
    out%ok=.true.; out%message='Panjer recursion'
  end function panjer_recursion

  function panjer_poisson(severity,lambda,tolerance,max_terms,scale) result(out)
    real(dp), intent(in) :: severity(:),lambda
    real(dp), intent(in), optional :: tolerance,scale
    integer, intent(in), optional :: max_terms
    type(aggregate_distribution) :: out
    real(dp) :: fs0
    fs0=exp(lambda*(severity(1)-1.0_dp))
    out=panjer_recursion(severity,0.0_dp,lambda,fs0,tolerance=tolerance, &
      max_terms=max_terms,scale=scale)
  end function panjer_poisson

  function panjer_binomial(severity,n,prob,tolerance,max_terms,scale) result(out)
    real(dp), intent(in) :: severity(:),prob
    integer, intent(in) :: n
    real(dp), intent(in), optional :: tolerance,scale
    integer, intent(in), optional :: max_terms
    type(aggregate_distribution) :: out
    real(dp) :: a,b,fs0
    a=-prob/(1.0_dp-prob); b=real(n+1,dp)*prob/(1.0_dp-prob)
    fs0=(1.0_dp+prob*(severity(1)-1.0_dp))**n
    out=panjer_recursion(severity,a,b,fs0,tolerance=tolerance, &
      max_terms=max_terms,scale=scale)
  end function panjer_binomial

  function panjer_nbinomial(severity,size,prob,tolerance,max_terms,scale) result(out)
    real(dp), intent(in) :: severity(:),size,prob
    real(dp), intent(in), optional :: tolerance,scale
    integer, intent(in), optional :: max_terms
    type(aggregate_distribution) :: out
    real(dp) :: a,b,fs0
    a=1.0_dp-prob; b=(size-1.0_dp)*a
    fs0=(prob/(1.0_dp-a*severity(1)))**size
    out=panjer_recursion(severity,a,b,fs0,tolerance=tolerance, &
      max_terms=max_terms,scale=scale)
  end function panjer_nbinomial

  function compound_simulation(nsim,frequency_rng,severity_rng,seed,step,max_value) result(out)
    integer, intent(in) :: nsim
    procedure(integer_rng_callback) :: frequency_rng
    procedure(real_rng_callback) :: severity_rng
    integer(i8), intent(in), optional :: seed
    real(dp), intent(in), optional :: step,max_value
    type(aggregate_distribution) :: out
    real(dp), allocatable :: values(:)
    real(dp) :: h,maxv,s
    integer :: i,j,n,bin,nbin
    if(present(seed)) call seed_rng(seed)
    h=1.0_dp; if(present(step)) h=step
    allocate(values(nsim))
    do i=1,nsim
      n=frequency_rng(); s=0.0_dp
      do j=1,max(0,n); s=s+severity_rng(); end do
      values(i)=s
    end do
    maxv=maxval(values); if(present(max_value)) maxv=max_value
    nbin=max(1,ceiling(maxv/h)+1)
    allocate(out%probability(nbin),out%support(nbin),out%cdf(nbin))
    out%probability=0.0_dp
    do i=1,nsim
      bin=min(nbin,max(1,nint(values(i)/h)+1))
      out%probability(bin)=out%probability(bin)+1.0_dp/real(nsim,dp)
    end do
    out%support=[(real(i-1,dp)*h,i=1,nbin)]
    out%cdf(1)=out%probability(1)
    do i=2,nbin; out%cdf(i)=out%cdf(i-1)+out%probability(i); end do
    out%ok=.true.; out%message='compound simulation'
  end function compound_simulation

  pure function aggregate_var(out,confidence) result(x)
    type(aggregate_distribution), intent(in) :: out
    real(dp), intent(in) :: confidence
    real(dp) :: x
    x=out%quantile(confidence)
  end function aggregate_var

  pure function aggregate_cte(out,confidence) result(x)
    type(aggregate_distribution), intent(in) :: out
    real(dp), intent(in) :: confidence
    real(dp) :: x,varv,den
    logical, allocatable :: mask(:)
    varv=out%quantile(confidence)
    allocate(mask(size(out%support))); mask=out%support>varv
    den=sum(out%probability,mask=mask)
    if(den<=0.0_dp) then
      x=varv
    else
      x=sum(out%support*out%probability,mask=mask)/den
    end if
  end function aggregate_cte

  pure function normal_aggregate_cdf(x,mean,variance) result(p)
    real(dp), intent(in) :: x,mean,variance
    real(dp) :: p
    if(variance<=0.0_dp) then
      p=merge(1.0_dp,0.0_dp,x>=mean)
    else
      p=normal_cdf((x-mean)/sqrt(variance))
    end if
  end function normal_aggregate_cdf

  pure function normal_power_cdf(x,mean,variance,skewness) result(p)
    real(dp), intent(in) :: x,mean,variance,skewness
    real(dp) :: p,z,y
    if(variance<=0.0_dp) then
      p=merge(1.0_dp,0.0_dp,x>=mean); return
    end if
    y=(x-mean)/sqrt(variance)
    if(abs(skewness)<1.0e-12_dp) then
      p=normal_cdf(y); return
    end if
    z=(-3.0_dp+sqrt(max(0.0_dp,9.0_dp+6.0_dp*skewness* &
      (y+skewness/6.0_dp))))/skewness
    p=normal_cdf(z)
  end function normal_power_cdf

end module actuar_aggregate
