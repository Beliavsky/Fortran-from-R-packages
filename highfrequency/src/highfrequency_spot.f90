! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_spot
  use highfrequency_kinds, only: dp, pi
  use highfrequency_stats, only: median3
  implicit none
  private
  public :: spot_volatility, spot_drift, drift_burst_statistic
  public :: preaverage_returns

contains

  pure real(dp) function kernel_value(kind,x) result(value)
    character(len=*),intent(in)::kind
    real(dp),intent(in)::x
    real(dp)::z
    z=abs(x)
    select case(trim(adjustl(kind)))
    case('gaussian','Gaussian')
      value=exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
    case('epanechnikov','Epanechnikov')
      value=0.75_dp*max(0.0_dp,1.0_dp-x*x)
    case('uniform','Uniform')
      value=merge(0.5_dp,0.0_dp,z<=1.0_dp)
    case('quadratic','Quadratic')
      value=0.9375_dp*max(0.0_dp,(1.0_dp-x*x)**2)
    case default
      value=exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
    end select
  end function kernel_value

  function spot_volatility(times,returns,evaluation_times,bandwidth,kernel) result(volatility)
    real(dp),intent(in)::times(:),returns(:),evaluation_times(:),bandwidth
    character(len=*),intent(in),optional::kernel
    real(dp)::volatility(size(evaluation_times))
    character(len=32)::kind
    real(dp)::w,denom
    integer::i,j,n
    n=min(size(times),size(returns))
    kind='gaussian';if(present(kernel))kind=kernel
    volatility=0.0_dp
    if(n==0 .or. bandwidth<=0.0_dp)return
    do j=1,size(evaluation_times)
      denom=0.0_dp
      do i=1,n
        w=kernel_value(kind,(times(i)-evaluation_times(j))/bandwidth)
        volatility(j)=volatility(j)+w*returns(i)*returns(i)
        denom=denom+w
      end do
      if(denom>0.0_dp)volatility(j)=volatility(j)/denom
    end do
  end function spot_volatility

  function spot_drift(times,returns,evaluation_times,bandwidth,method,kernel) result(drift)
    real(dp),intent(in)::times(:),returns(:),evaluation_times(:),bandwidth
    character(len=*),intent(in),optional::method,kernel
    real(dp)::drift(size(evaluation_times))
    character(len=32)::meth,kind
    real(dp),allocatable::values(:),weights(:)
    real(dp)::w,denom
    integer::i,j,n,m
    n=min(size(times),size(returns))
    meth='mean';if(present(method))meth=method
    kind='gaussian';if(present(kernel))kind=kernel
    drift=0.0_dp
    if(n==0 .or. bandwidth<=0.0_dp)return
    allocate(values(n),weights(n))
    do j=1,size(evaluation_times)
      m=0
      denom=0.0_dp
      do i=1,n
        w=kernel_value(kind,(times(i)-evaluation_times(j))/bandwidth)
        if(w>0.0_dp)then
          m=m+1
          values(m)=returns(i)
          weights(m)=w
          denom=denom+w
        end if
      end do
      if(m==0)cycle
      if(trim(adjustl(meth))=='median')then
        call sort_values(values(:m))
        if(mod(m,2)==1)then
          drift(j)=values((m+1)/2)
        else
          drift(j)=0.5_dp*(values(m/2)+values(m/2+1))
        end if
      else
        drift(j)=sum(values(:m)*weights(:m))/max(denom,tiny(1.0_dp))
      end if
    end do
  contains
    subroutine sort_values(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::ii,jj
      do ii=2,size(x)
        key=x(ii);jj=ii-1
        do while(jj>=1)
          if(x(jj)<=key)exit
          x(jj+1)=x(jj);jj=jj-1
        end do
        x(jj+1)=key
      end do
    end subroutine sort_values
  end function spot_drift

  function drift_burst_statistic(times,returns,evaluation_times,bandwidth,kernel) result(statistic)
    real(dp),intent(in)::times(:),returns(:),evaluation_times(:),bandwidth
    character(len=*),intent(in),optional::kernel
    real(dp)::statistic(size(evaluation_times))
    real(dp),allocatable::drift(:),variance(:)
    allocate(drift(size(evaluation_times)),variance(size(evaluation_times)))
    drift=spot_drift(times,returns,evaluation_times,bandwidth,'mean',kernel)
    variance=spot_volatility(times,returns,evaluation_times,bandwidth,kernel)
    statistic=drift/sqrt(max(variance,tiny(1.0_dp)))
  end function drift_burst_statistic

  function preaverage_returns(returns,window) result(preaveraged)
    real(dp),intent(in)::returns(:)
    integer,intent(in)::window
    real(dp),allocatable::preaveraged(:)
    integer::i,j,n,m
    real(dp)::weight
    n=size(returns)
    if(window<2 .or. n<window)then
      allocate(preaveraged(0))
      return
    end if
    m=n-window+1
    allocate(preaveraged(m))
    preaveraged=0.0_dp
    do i=1,m
      do j=1,window-1
        weight=min(real(j,dp)/real(window,dp),1.0_dp-real(j,dp)/real(window,dp))
        preaveraged(i)=preaveraged(i)+weight*returns(i+j-1)
      end do
    end do
  end function preaverage_returns

end module highfrequency_spot
