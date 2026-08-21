! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_statistics
  use tsa_kinds, only : dp
  use tsa_types, only : spectrum_result, runs_result
  use tsa_utils, only : mean_value, variance_n, sd_n
  implicit none
  private
  public :: skewness, kurtosis, autocorrelation, autocovariance, partial_autocorrelation
  public :: cross_correlation, periodogram, harmonic_matrix, season_index, runs

contains

  pure real(dp) function skewness(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m,s
    if(size(x)==0) then
    value=0.0_dp
    return
    end if
    m=mean_value(x)
    s=sd_n(x)
    if(s<=tiny(1.0_dp)) then
    value=0.0_dp
    else
    value=sum((x-m)**3)/(real(size(x),dp)*s**3)
    end if
  end function skewness

  pure real(dp) function kurtosis(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m,v
    if(size(x)==0) then
    value=0.0_dp
    return
    end if
    m=mean_value(x)
    v=variance_n(x)
    if(v<=tiny(1.0_dp)) then
    value=0.0_dp
    else
    value=sum((x-m)**4)/(real(size(x),dp)*v*v)-3.0_dp
    end if
  end function kurtosis

  subroutine autocovariance(x,lag_max,cov,demean)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag_max
    real(dp), intent(out) :: cov(0:)
    logical, intent(in), optional :: demean
    logical :: dm
    real(dp) :: m
    integer :: k,n
    dm=.true.
    if(present(demean)) dm=demean
    n=size(x)
    m=merge(mean_value(x),0.0_dp,dm)
    cov=0.0_dp
    do k=0,min(lag_max,ubound(cov,1))
      cov(k)=sum((x(1:n-k)-m)*(x(1+k:n)-m))/real(n,dp)
    end do
  end subroutine autocovariance

  subroutine autocorrelation(x,lag_max,acf,demean,drop_lag_zero)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag_max
    real(dp), allocatable, intent(out) :: acf(:)
    logical, intent(in), optional :: demean,drop_lag_zero
    real(dp), allocatable :: cv(:)
    logical :: drop
    integer :: kmax
    kmax=min(lag_max,max(0,size(x)-1))
    allocate(cv(0:kmax))
    call autocovariance(x,kmax,cv,demean)
    if(cv(0)<=tiny(1.0_dp)) then
      cv=0.0_dp
      cv(0)=1.0_dp
    else
      cv=cv/cv(0)
    end if
    drop=.true.
    if(present(drop_lag_zero)) drop=drop_lag_zero
    if(drop) then
      allocate(acf(kmax))
      if(kmax>0) acf=cv(1:)
    else
      allocate(acf(kmax+1))
      acf=cv
    end if
  end subroutine autocorrelation

  subroutine partial_autocorrelation(x,lag_max,pacf)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: lag_max
    real(dp), allocatable, intent(out) :: pacf(:)
    real(dp), allocatable :: ac(:),phi(:,:),v(:)
    integer :: k,j,m
    m=min(lag_max,max(0,size(x)-2))
    allocate(pacf(m))
    if(m==0)return
    call autocorrelation(x,m,ac,drop_lag_zero=.false.)
    allocate(phi(m,m),v(m))
    phi=0.0_dp
    phi(1,1)=ac(2)
    pacf(1)=phi(1,1)
    v(1)=1.0_dp-phi(1,1)**2
    do k=2,m
      phi(k,k)=ac(k+1)
      do j=1,k-1
        phi(k,k)=phi(k,k)-phi(k-1,j)*ac(k-j+1)
      end do
      phi(k,k)=phi(k,k)/max(v(k-1),tiny(1.0_dp))
      do j=1,k-1
        phi(k,j)=phi(k-1,j)-phi(k,k)*phi(k-1,k-j)
      end do
      v(k)=v(k-1)*(1.0_dp-phi(k,k)**2)
      pacf(k)=phi(k,k)
    end do
  end subroutine partial_autocorrelation

  subroutine cross_correlation(x,y,lag_max,lags,ccf)
    real(dp), intent(in) :: x(:),y(:)
    integer, intent(in) :: lag_max
    integer, allocatable, intent(out) :: lags(:)
    real(dp), allocatable, intent(out) :: ccf(:)
    integer :: k,n,i,nuse
    real(dp) :: mx,my,sx,sy
    n=min(size(x),size(y))
    mx=mean_value(x(:n))
    my=mean_value(y(:n))
    sx=sd_n(x(:n))
    sy=sd_n(y(:n))
    allocate(lags(2*lag_max+1),ccf(2*lag_max+1))
    i=0
    do k=-lag_max,lag_max
      i=i+1
      lags(i)=k
      if(k>=0) then
        nuse=n-k
        ccf(i)=sum((x(1:nuse)-mx)*(y(1+k:n)-my))
      else
        nuse=n+k
        ccf(i)=sum((x(1-k:n)-mx)*(y(1:nuse)-my))
      end if
      if(nuse>0 .and. sx>0.0_dp .and. sy>0.0_dp) then
        ccf(i)=ccf(i)/(real(n,dp)*sx*sy)
      else
        ccf(i)=0.0_dp
      end if
    end do
  end subroutine cross_correlation

  function periodogram(x,demean,detrend) result(res)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: demean,detrend
    type(spectrum_result) :: res
    real(dp), allocatable :: y(:)
    logical :: dm,dt
    integer :: n,k,j,m
    real(dp) :: pi,ang,re,im,tmean,den,slope,intercept
    n=size(x)
    dm=.true.
    dt=.false.
    if(present(demean))dm=demean
    if(present(detrend))dt=detrend
    if(n<2) then
    res%status=1
    allocate(res%frequency(0),res%spectrum(0))
    return
    end if
    allocate(y(n))
    y=x
    if(dt) then
      tmean=0.5_dp*real(n+1,dp)
      den=sum([( (real(j,dp)-tmean)**2,j=1,n )])
      slope=sum([( (real(j,dp)-tmean)*(y(j)-mean_value(y)),j=1,n )])/max(den,tiny(1.0_dp))
      intercept=mean_value(y)-slope*tmean
      do j=1,n
      y(j)=y(j)-(intercept+slope*real(j,dp))
      end do
    else if(dm) then
      y=y-mean_value(y)
    end if
    m=n/2
    allocate(res%frequency(m),res%spectrum(m))
    pi=acos(-1.0_dp)
    do k=1,m
      re=0.0_dp
      im=0.0_dp
      do j=1,n
        ang=2.0_dp*pi*real(k*j,dp)/real(n,dp)
        re=re+y(j)*cos(ang)
        im=im-y(j)*sin(ang)
      end do
      res%frequency(k)=real(k,dp)/real(n,dp)
      res%spectrum(k)=(re*re+im*im)/real(n,dp)
      ! TSA periodogram doubles the one-sided spectrum except at Nyquist.
      if(.not.(mod(n,2)==0 .and. k==m)) res%spectrum(k)=2.0_dp*res%spectrum(k)
    end do
  end function periodogram

  function harmonic_matrix(n,frequency,m,start_cycle) result(h)
    integer, intent(in) :: n,frequency,m
    real(dp), intent(in), optional :: start_cycle
    real(dp), allocatable :: h(:,:)
    real(dp) :: t0,t,pi
    integer :: i,j,ncol,col
    if(2*m>frequency .or. m<1) then
    allocate(h(0,0))
    return
    end if
    ncol=2*m
    if(2*m==frequency)ncol=ncol-1
    allocate(h(n,ncol))
    pi=acos(-1.0_dp)
    t0=0.0_dp
    if(present(start_cycle))t0=start_cycle
    do i=1,n
      t=t0+real(i-1,dp)/real(frequency,dp)
      do j=1,m
        h(i,j)=cos(2.0_dp*pi*real(j,dp)*t)
        col=m+j
        if(col<=ncol)h(i,col)=sin(2.0_dp*pi*real(j,dp)*t)
      end do
    end do
  end function harmonic_matrix

  function season_index(n,frequency,start_season) result(s)
    integer, intent(in) :: n,frequency
    integer, intent(in), optional :: start_season
    integer, allocatable :: s(:)
    integer :: i,st
    st=1
    if(present(start_season))st=start_season
    allocate(s(n))
    do i=1,n
    s(i)=1+mod(st-1+i-1,frequency)
    end do
  end function season_index

  function runs(x,k) result(res)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: k
    type(runs_result) :: res
    integer, allocatable :: y(:)
    integer :: i,r,l2,j,rr
    real(dp) :: th,mu,pv
    real(dp), allocatable :: pdf(:),f(:),g1(:),g2(:)
    th=0.0_dp
    if(present(k))th=k
    res%threshold=th
    allocate(y(size(x)))
    y=merge(1,0,x<=th)
    res%n1=sum(y)
    res%n2=size(y)-res%n1
    res%expected_runs=1.0_dp+2.0_dp*real(res%n1*res%n2,dp)/real(max(1,res%n1+res%n2),dp)
    if(res%n1*res%n2==0) then
    res%observed_runs=1
    res%p_value=-1.0_dp
    return
    end if
    r=1
    do i=2,size(y)
    if(y(i)/=y(i-1))r=r+1
    end do
    res%observed_runs=r
    if(res%n1==res%n2) then
    l2=2*res%n1
    else
    l2=2*min(res%n1,res%n2)+1
    end if
    allocate(pdf(l2),f(max(1,l2/2)),g1(max(1,l2/2)),g2(max(1,l2/2)))
    pdf=0.0_dp
    f=0.0_dp
    g1=0.0_dp
    g2=0.0_dp
    f(1)=2.0_dp
    g1(1)=real(res%n1-1,dp)
    g2(1)=real(res%n2-1,dp)
    if(l2>=2)pdf(2)=f(1)
    if(l2>=3)pdf(3)=g1(1)+g2(1)
    if(l2>4) then
      do i=4,l2,2
        rr=(i-2)/2
        j=rr+1
        f(j)=real((res%n1-rr)*(res%n2-rr),dp)/real(rr*rr,dp)*f(rr)
        pdf(i)=f(j)
      end do
    end if
    if(l2>5) then
      do i=5,l2,2
        rr=(i-3)/2
        j=rr+1
        g1(j)=real((res%n1-rr-1)*(res%n2-rr),dp)/real((rr+1)*rr,dp)*g1(rr)
        g2(j)=real((res%n2-rr-1)*(res%n1-rr),dp)/real((rr+1)*rr,dp)*g2(rr)
        pdf(i)=g1(j)+g2(j)
      end do
    end if
    pdf=pdf/sum(pdf)
    mu=res%expected_runs
    if(real(r,dp)<=mu) then
    pv=sum(pdf,mask=[(i<=r,i=1,l2)])
    else
    pv=sum(pdf,mask=[(i>=r,i=1,l2)])
    end if
    if(pv>0.5_dp)pv=1.0_dp-pv
    res%p_value=min(1.0_dp,2.0_dp*pv)
  end function runs
end module tsa_statistics
