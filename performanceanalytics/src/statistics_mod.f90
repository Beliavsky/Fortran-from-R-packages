! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module statistics_mod
  use kinds_mod, only: dp
  implicit none
  private
  public :: mean_value, variance_value, sd_value, covariance_value, correlation_value
  public :: median_value, quantile_type7, centered_moment, skewness_value, kurtosis_value
  public :: mean_absolute_deviation, geometric_mean, standard_error, autocorrelation
  public :: lower_partial_moment, upper_partial_moment, percentile_rank
  public :: sort_real, solve_linear_system, ols_fit
contains
  pure real(dp) function mean_value(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = 0.0_dp
    else
      v = sum(x)/real(size(x),dp)
    end if
  end function mean_value

  pure real(dp) function variance_value(x, sample) result(v)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: sample
    logical :: use_sample
    real(dp) :: m, den
    use_sample = .true.
    if (present(sample)) use_sample = sample
    if (size(x) <= merge(1,0,use_sample)) then
      v = 0.0_dp
      return
    end if
    m = mean_value(x)
    den = real(size(x)-merge(1,0,use_sample),dp)
    v = sum((x-m)**2)/den
  end function variance_value

  pure real(dp) function sd_value(x, sample) result(v)
    real(dp), intent(in) :: x(:)
    logical, intent(in), optional :: sample
    if (present(sample)) then
      v = sqrt(max(0.0_dp,variance_value(x,sample)))
    else
      v = sqrt(max(0.0_dp,variance_value(x)))
    end if
  end function sd_value

  pure real(dp) function covariance_value(x,y,sample) result(v)
    real(dp), intent(in) :: x(:), y(:)
    logical, intent(in), optional :: sample
    logical :: use_sample
    real(dp) :: den
    integer :: n
    n = min(size(x),size(y))
    use_sample=.true.; if (present(sample)) use_sample=sample
    if (n <= merge(1,0,use_sample)) then
      v=0.0_dp; return
    end if
    den=real(n-merge(1,0,use_sample),dp)
    v=sum((x(:n)-mean_value(x(:n)))*(y(:n)-mean_value(y(:n))))/den
  end function covariance_value

  pure real(dp) function correlation_value(x,y) result(v)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: sx,sy
    sx=sd_value(x); sy=sd_value(y)
    if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
      v=0.0_dp
    else
      v=covariance_value(x,y)/(sx*sy)
    end if
  end function correlation_value

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i); j=i-1
      do while(j>=1)
        if (x(j)<=key) exit
        x(j+1)=x(j); j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  real(dp) function quantile_type7(x,p) result(q)
    real(dp), intent(in) :: x(:),p
    real(dp), allocatable :: z(:)
    real(dp) :: h,g
    integer :: j,n
    n=size(x)
    if(n==0) then; q=0.0_dp; return; end if
    allocate(z(n)); z=x; call sort_real(z)
    if(p<=0.0_dp) then; q=z(1); return; end if
    if(p>=1.0_dp) then; q=z(n); return; end if
    h=1.0_dp+real(n-1,dp)*p; j=int(floor(h)); g=h-real(j,dp)
    q=(1.0_dp-g)*z(j)+g*z(min(j+1,n))
  end function quantile_type7

  real(dp) function median_value(x) result(v)
    real(dp), intent(in) :: x(:)
    v=quantile_type7(x,0.5_dp)
  end function median_value

  pure real(dp) function centered_moment(x,power) result(v)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: power
    real(dp) :: m
    if(size(x)==0) then; v=0.0_dp; return; end if
    m=mean_value(x); v=sum((x-m)**power)/real(size(x),dp)
  end function centered_moment

  pure real(dp) function skewness_value(x,method) result(v)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: method
    integer :: n,meth
    real(dp) :: m2,m3,g1
    n=size(x); meth=3; if(present(method)) meth=method
    m2=centered_moment(x,2); m3=centered_moment(x,3)
    if(m2<=tiny(1.0_dp)) then; v=0.0_dp; return; end if
    g1=m3/m2**1.5_dp
    select case(meth)
    case(1); v=g1
    case(2)
      if(n>2) then; v=sqrt(real(n*(n-1),dp))/real(n-2,dp)*g1; else; v=0.0_dp; end if
    case default
      if(n>1) then; v=g1*(real(n-1,dp)/real(n,dp))**1.5_dp; else; v=0.0_dp; end if
    end select
  end function skewness_value

  pure real(dp) function kurtosis_value(x,method,excess) result(v)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: method
    logical, intent(in), optional :: excess
    integer :: n,meth
    logical :: ex
    real(dp) :: m2,m4,g2
    n=size(x); meth=3; ex=.true.
    if(present(method)) meth=method
    if(present(excess)) ex=excess
    m2=centered_moment(x,2); m4=centered_moment(x,4)
    if(m2<=tiny(1.0_dp)) then; v=0.0_dp; return; end if
    g2=m4/(m2*m2)-3.0_dp
    select case(meth)
    case(1); v=g2
    case(2)
      if(n>3) then
        v=real(n-1,dp)/real((n-2)*(n-3),dp)*((real(n+1,dp))*g2+6.0_dp)
      else
        v=0.0_dp
      end if
    case default
      if(n>1) then
        v=(real(n-1,dp)/real(n,dp))**2*(g2+3.0_dp)-3.0_dp
      else
        v=0.0_dp
      end if
    end select
    if(.not.ex) v=v+3.0_dp
  end function kurtosis_value

  pure real(dp) function mean_absolute_deviation(x,center) result(v)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: center
    real(dp) :: c
    if(size(x)==0) then; v=0.0_dp; return; end if
    c=mean_value(x); if(present(center)) c=center
    v=sum(abs(x-c))/real(size(x),dp)
  end function mean_absolute_deviation

  pure real(dp) function geometric_mean(x) result(v)
    real(dp), intent(in) :: x(:)
    if(size(x)==0) then; v=0.0_dp; return; end if
    if(any(1.0_dp+x<=0.0_dp)) then
      v=-1.0_dp
    else
      v=exp(sum(log(1.0_dp+x))/real(size(x),dp))-1.0_dp
    end if
  end function geometric_mean

  pure real(dp) function standard_error(x) result(v)
    real(dp), intent(in) :: x(:)
    if(size(x)==0) then; v=0.0_dp; else; v=sd_value(x)/sqrt(real(size(x),dp)); end if
  end function standard_error

  pure real(dp) function autocorrelation(x,lag) result(v)
    real(dp), intent(in) :: x(:)
    integer,intent(in)::lag
    integer::n
    n=size(x)
    if(lag<0 .or. lag>=n) then; v=0.0_dp; return; end if
    if(lag==0) then; v=1.0_dp; else; v=correlation_value(x(1:n-lag),x(1+lag:n)); end if
  end function autocorrelation

  pure real(dp) function lower_partial_moment(x,threshold,order) result(v)
    real(dp),intent(in)::x(:),threshold
    integer,intent(in)::order
    if(size(x)==0) then; v=0.0_dp; else; v=sum(max(threshold-x,0.0_dp)**order)/real(size(x),dp); end if
  end function lower_partial_moment

  pure real(dp) function upper_partial_moment(x,threshold,order) result(v)
    real(dp),intent(in)::x(:),threshold
    integer,intent(in)::order
    if(size(x)==0) then; v=0.0_dp; else; v=sum(max(x-threshold,0.0_dp)**order)/real(size(x),dp); end if
  end function upper_partial_moment

  real(dp) function percentile_rank(x,value) result(v)
    real(dp),intent(in)::x(:),value
    if(size(x)==0) then; v=0.0_dp; else; v=real(count(x<=value),dp)/real(size(x),dp); end if
  end function percentile_rank

  subroutine solve_linear_system(a,b,x,ok)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),intent(out)::x(:)
    logical,intent(out)::ok
    real(dp),allocatable::aa(:,:),bb(:)
    real(dp)::fac,piv,tmp
    integer::n,i,j,k,imax
    n=size(b); ok=.false.; x=0.0_dp
    if(size(a,1)/=n .or. size(a,2)/=n .or. size(x)/=n) return
    allocate(aa(n,n),bb(n)); aa=a; bb=b
    do k=1,n-1
      imax=k
      do i=k+1,n
        if(abs(aa(i,k))>abs(aa(imax,k))) imax=i
      end do
      if(abs(aa(imax,k))<=1.0e-14_dp) return
      if(imax/=k) then
        do j=k,n; tmp=aa(k,j); aa(k,j)=aa(imax,j); aa(imax,j)=tmp; end do
        tmp=bb(k); bb(k)=bb(imax); bb(imax)=tmp
      end if
      do i=k+1,n
        fac=aa(i,k)/aa(k,k)
        aa(i,k:n)=aa(i,k:n)-fac*aa(k,k:n); bb(i)=bb(i)-fac*bb(k)
      end do
    end do
    if(abs(aa(n,n))<=1.0e-14_dp) return
    do i=n,1,-1
      piv=bb(i)-sum(aa(i,i+1:n)*x(i+1:n)); x(i)=piv/aa(i,i)
    end do
    ok=.true.
  end subroutine solve_linear_system

  subroutine ols_fit(x,y,beta,residuals,r2,ok)
    real(dp),intent(in)::x(:,:),y(:)
    real(dp),intent(out)::beta(:),residuals(:),r2
    logical,intent(out)::ok
    real(dp),allocatable::xtx(:,:),xty(:)
    real(dp)::sst
    integer::n,p
    n=size(y); p=size(x,2); ok=.false.; r2=0.0_dp
    if(size(x,1)/=n .or. size(beta)/=p .or. size(residuals)/=n) return
    allocate(xtx(p,p),xty(p)); xtx=matmul(transpose(x),x); xty=matmul(transpose(x),y)
    call solve_linear_system(xtx,xty,beta,ok); if(.not.ok) return
    residuals=y-matmul(x,beta); sst=sum((y-mean_value(y))**2)
    if(sst>tiny(1.0_dp)) r2=1.0_dp-sum(residuals**2)/sst
  end subroutine ols_fit
end module statistics_mod
