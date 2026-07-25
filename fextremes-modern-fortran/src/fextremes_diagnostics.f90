! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_diagnostics
  use fextremes_kinds, only: dp
  use fextremes_stats, only: sort_real, mean_value, variance_value, normal_quantile, autocorrelation, &
    finite_mean_sd, nan_value
  use fextremes_distributions, only: gpd_quantile
  implicit none
  private
  public :: curve_result, record_result, exceedance_acf_result, tail_index_result
  public :: empirical_survival, pareto_qq, mean_excess_curve, mean_residual_life
  public :: records_development, subsample_record_counts, max_sum_ratios, slln_path, lil_path
  public :: exceedance_acf, pickands_estimator, hill_estimator, dehaan_estimator
  public :: normal_mean_excess

  type :: curve_result
    real(dp), allocatable :: x(:), y(:), lower(:), upper(:)
  end type curve_result
  type :: record_result
    integer, allocatable :: number(:), trial(:)
    real(dp), allocatable :: record(:), expected(:), se(:)
  end type record_result
  type :: exceedance_acf_result
    real(dp), allocatable :: heights(:), distances(:), height_acf(:), distance_acf(:)
  end type exceedance_acf_result
  type :: tail_index_result
    real(dp), allocatable :: estimates(:)
    real(dp) :: mean=0.0_dp, sd=0.0_dp
  end type tail_index_result
contains
  subroutine empirical_survival(data,result)
    real(dp),intent(in)::data(:); type(curve_result),intent(out)::result
    integer::i,n
    n=size(data); allocate(result%x(n),result%y(n)); result%x=data; call sort_real(result%x)
    do i=1,n; result%y(i)=1.0_dp-(real(i,dp)-0.5_dp)/real(n,dp); end do
  end subroutine empirical_survival

  subroutine pareto_qq(data,xi,result)
    real(dp),intent(in)::data(:),xi; type(curve_result),intent(out)::result
    integer::i,n; real(dp)::p
    n=size(data); allocate(result%x(n),result%y(n)); result%x=data; call sort_real(result%x)
    do i=1,n
      p=(real(i,dp)-0.5_dp)/real(n,dp)
      if(abs(xi)<1.0e-10_dp) then; result%y(i)=-log(1.0_dp-p)
      else; result%y(i)=gpd_quantile(p,xi,0.0_dp,1.0_dp); end if
    end do
  end subroutine pareto_qq

  subroutine mean_excess_curve(data,result)
    real(dp),intent(in)::data(:); type(curve_result),intent(out)::result
    real(dp),allocatable::s(:)
    integer::i,j,n,m
    n=size(data); allocate(s(n)); s=data; call sort_real(s)
    m=1
    do i=2,n-1
      if(abs(s(i)-s(i-1))>0.0_dp) m=m+1
    end do
    allocate(result%x(m),result%y(m))
    j=1
    result%x(j)=s(1)
    result%y(j)=mean_value(s(2:n)-s(1))
    do i=2,n-1
      if(abs(s(i)-s(i-1))>0.0_dp) then
        j=j+1
        result%x(j)=s(i)
        result%y(j)=mean_value(s(i+1:n)-s(i))
      end if
    end do
  end subroutine mean_excess_curve

  subroutine mean_residual_life(data,umin,umax,npoints,confidence,result)
    real(dp),intent(in)::data(:),umin,umax,confidence
    integer,intent(in)::npoints; type(curve_result),intent(out)::result
    real(dp),allocatable::excess(:); real(dp)::u,z,sd
    integer::i,j,n
    allocate(result%x(npoints),result%y(npoints),result%lower(npoints),result%upper(npoints))
    z=normal_quantile(0.5_dp*(1.0_dp+confidence))
    do i=1,npoints
      if(npoints==1) then; u=umin; else; u=umin+(umax-umin)*real(i-1,dp)/real(npoints-1,dp); end if
      n=count(data>=u); result%x(i)=u
      if(n<2) then; result%y(i)=nan_value(); result%lower(i)=nan_value(); result%upper(i)=nan_value(); cycle; end if
      allocate(excess(n)); j=0
      do concurrent (j=1:n); excess(j)=0.0_dp; end do
      j=0
      block
        integer :: k
        do k=1,size(data); if(data(k)>=u) then; j=j+1; excess(j)=data(k)-u; end if; end do
      end block
      result%y(i)=mean_value(excess); sd=sqrt(variance_value(excess))
      result%lower(i)=result%y(i)-z*sd/sqrt(real(n,dp)); result%upper(i)=result%y(i)+z*sd/sqrt(real(n,dp))
      deallocate(excess)
    end do
  end subroutine mean_residual_life

  subroutine records_development(data,result)
    real(dp),intent(in)::data(:); type(record_result),intent(out)::result
    real(dp),allocatable::cum(:); integer::i,j,m
    allocate(cum(size(data))); cum(1)=data(1); do i=2,size(data); cum(i)=max(cum(i-1),data(i)); end do
    m=1; do i=2,size(data); if(cum(i)>cum(i-1)) m=m+1; end do
    allocate(result%number(m),result%trial(m),result%record(m),result%expected(m),result%se(m))
    m=1
    result%number(m)=m
    result%trial(m)=1
    result%record(m)=cum(1)
    result%expected(m)=1.0_dp
    result%se(m)=0.0_dp
    do i=2,size(data)
      if(cum(i)>cum(i-1)) then
        m=m+1
        result%number(m)=m
        result%trial(m)=i
        result%record(m)=cum(i)
        result%expected(m)=sum(1.0_dp/real([(j,j=1,i)],dp))
        result%se(m)=sqrt(max(0.0_dp,result%expected(m)- &
          sum(1.0_dp/real([(j*j,j=1,i)],dp))))
      end if
    end do
  end subroutine records_development

  subroutine subsample_record_counts(data,nsub,counts)
    real(dp),intent(in)::data(:); integer,intent(in)::nsub; integer,allocatable,intent(out)::counts(:)
    integer::len,i,j,lo,hi,c; real(dp)::rmax
    len=size(data)/nsub; allocate(counts(nsub)); counts=0
    do i=1,nsub
      lo=(i-1)*len+1; hi=i*len; if(lo>hi) cycle
      c=1; rmax=data(lo)
      do j=lo+1,hi; if(data(j)>rmax) then; c=c+1; rmax=data(j); end if; end do
      counts(i)=c
    end do
  end subroutine subsample_record_counts

  subroutine max_sum_ratios(data,powers,ratios)
    real(dp),intent(in)::data(:),powers(:); real(dp),allocatable,intent(out)::ratios(:,:)
    integer::i,j; real(dp)::mx,sm,v
    allocate(ratios(size(data),size(powers)))
    do j=1,size(powers)
      mx=0.0_dp; sm=0.0_dp
      do i=1,size(data); v=abs(data(i))**powers(j); mx=max(mx,v); sm=sm+v; ratios(i,j)=mx/sm; end do
    end do
  end subroutine max_sum_ratios

  subroutine slln_path(data,path)
    real(dp),intent(in)::data(:); real(dp),allocatable,intent(out)::path(:)
    integer::i; allocate(path(size(data))); path(1)=data(1)
    do i=2,size(data); path(i)=(path(i-1)*real(i-1,dp)+data(i))/real(i,dp); end do
  end subroutine slln_path

  subroutine lil_path(data,path)
    real(dp),intent(in)::data(:); real(dp),allocatable,intent(out)::path(:)
    real(dp)::m,sd,cs; integer::i
    allocate(path(max(0,size(data)-2))); if(size(data)<3) return
    m=mean_value(data); sd=sqrt(variance_value(data)); cs=sum(data(1:2))
    do i=3,size(data); cs=cs+data(i); path(i-2)=(cs-m*real(i,dp))/(sqrt(2.0_dp*real(i,dp)*log(log(real(i,dp))))*sd); end do
  end subroutine lil_path

  subroutine exceedance_acf(data,threshold,max_lag,result)
    real(dp),intent(in)::data(:),threshold; integer,intent(in)::max_lag
    type(exceedance_acf_result),intent(out)::result
    integer,allocatable::idx(:); integer::i,j,n
    n=count(data>threshold); allocate(result%heights(n),idx(n)); j=0
    do i=1,size(data); if(data(i)>threshold) then; j=j+1; result%heights(j)=data(i)-threshold; idx(j)=i; end if; end do
    allocate(result%distances(max(0,n-1)),result%height_acf(0:max_lag),result%distance_acf(0:max_lag))
    if(n>1) result%distances=real(idx(2:)-idx(:n-1),dp)
    call autocorrelation(result%heights,max_lag,result%height_acf)
    if(n>2) then; call autocorrelation(result%distances,max_lag,result%distance_acf); else; result%distance_acf=nan_value(); end if
  end subroutine exceedance_acf

  subroutine pickands_estimator(values,fraction,result)
    real(dp),intent(in)::values(:),fraction; type(tail_index_result),intent(out)::result
    real(dp),allocatable::s(:); integer::n,k,i,nvalid
    n=max(4,min(size(values),floor(fraction*real(size(values),dp))))
    allocate(s(size(values)))
    s=abs(values)
    call sort_real(s,ascending=.false.)
    k=n/4; allocate(result%estimates(k))
    do i=1,k
      if((s(i)-s(2*i))>0.0_dp .and. (s(2*i)-s(4*i))>0.0_dp) then
        result%estimates(i)=log((s(i)-s(2*i))/(s(2*i)-s(4*i)))/log(2.0_dp)
      else; result%estimates(i)=nan_value(); end if
    end do
    call finite_mean_sd(result%estimates,result%mean,result%sd,nvalid)
  end subroutine pickands_estimator

  subroutine hill_estimator(values,fraction,result)
    real(dp),intent(in)::values(:),fraction; type(tail_index_result),intent(out)::result
    real(dp),allocatable::s(:); real(dp)::cs; integer::n,i,nvalid
    n=max(2,min(size(values),floor(fraction*real(size(values),dp))))
    allocate(s(size(values)))
    s=abs(values)
    call sort_real(s,ascending=.false.)
    allocate(result%estimates(n-1)); cs=log(s(1))
    do i=2,n; cs=cs+log(s(i)); result%estimates(i-1)=cs/real(i,dp)-log(s(i)); end do
    call finite_mean_sd(result%estimates,result%mean,result%sd,nvalid)
  end subroutine hill_estimator

  subroutine dehaan_estimator(values,fraction,result)
    real(dp),intent(in)::values(:),fraction; type(tail_index_result),intent(out)::result
    real(dp),allocatable::s(:); real(dp)::cs,cs2,lx,m1,m2; integer::n,i,nvalid
    n=max(2,min(size(values),floor(fraction*real(size(values),dp))))
    allocate(s(size(values)))
    s=abs(values)
    call sort_real(s,ascending=.false.)
    allocate(result%estimates(n-1)); cs=0.0_dp; cs2=0.0_dp
    do i=1,n-1
      lx=log(s(i)); cs=cs+lx; cs2=cs2+lx*lx
      m1=cs/real(i,dp)-log(s(i+1)); m2=cs2/real(i,dp)-2.0_dp*cs*log(s(i+1))/real(i,dp)+log(s(i+1))**2
      if(m2>0.0_dp .and. abs(m1*m1/m2-1.0_dp)>1.0e-12_dp) then; result%estimates(i)=1.0_dp+m1+0.5_dp/(m1*m1/m2-1.0_dp)
      else; result%estimates(i)=nan_value(); end if
    end do
    call finite_mean_sd(result%estimates,result%mean,result%sd,nvalid)
  end subroutine dehaan_estimator

  pure real(dp) function normal_mean_excess(threshold,mu,sigma) result(me)
    real(dp),intent(in)::threshold,mu,sigma; real(dp)::z,surv
    if(sigma<=0.0_dp) then; me=nan_value(); return; end if
    z=(threshold-mu)/sigma; surv=0.5_dp*erfc(z/sqrt(2.0_dp))
    if(surv<=tiny(1.0_dp)) then; me=0.0_dp; else; me=mu-threshold+sigma*exp(-0.5_dp*z*z)/sqrt(2.0_dp*acos(-1.0_dp))/surv; end if
  end function normal_mean_excess
end module fextremes_diagnostics
