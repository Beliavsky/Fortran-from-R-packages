! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_utils
  use tsdyn_kinds, only: dp, pi_dp, n_deterministic, include_none, include_const, include_trend, include_both
  use tsdyn_linalg, only: cholesky_lower
  implicit none
  private
  public :: seed_random, random_normal, random_normal_vector, random_normal_matrix
  public :: build_deterministic, lag_embed_univariate, lag_embed_multivariate
  public :: differences, quantile_linear, sorted_unique, argsort_real
  public :: sigmoid, normal_cdf, normal_pdf, sample_with_replacement
  public :: mean_value, variance_value, percentile_grid
contains
  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n,i
    integer, allocatable :: put(:)
    call random_seed(size=n); allocate(put(n))
    do i=1,n
      put(i)=mod(abs(seed)+104729*i,2147483646)+1
    end do
    call random_seed(put=put)
  end subroutine seed_random

  real(dp) function random_normal() result(z)
    real(dp) :: u1,u2
    call random_number(u1); call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi_dp*u2)
  end function random_normal

  subroutine random_normal_vector(z)
    real(dp), intent(out) :: z(:)
    integer::i
    do i=1,size(z); z(i)=random_normal(); end do
  end subroutine random_normal_vector

  subroutine random_normal_matrix(z, sigma, info)
    real(dp), intent(out) :: z(:,:)
    real(dp), intent(in), optional :: sigma(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: l(:,:),u(:)
    integer::i,j,istat
    if(present(sigma))then
      call cholesky_lower(sigma,l,istat)
      if(present(info))info=istat
      if(istat/=0)then; z=0.0_dp; return; end if
      allocate(u(size(z,2)))
      do i=1,size(z,1)
        call random_normal_vector(u)
        z(i,:)=matmul(l,u)
      end do
    else
      if(present(info))info=0
      do j=1,size(z,2); call random_normal_vector(z(:,j)); end do
    end if
  end subroutine random_normal_matrix

  subroutine build_deterministic(n, include, dmat, trend_start)
    integer,intent(in)::n,include
    real(dp),allocatable,intent(out)::dmat(:,:)
    real(dp),intent(in),optional::trend_start
    integer::i,nd
    real(dp)::s
    nd=n_deterministic(include); s=1.0_dp; if(present(trend_start))s=trend_start
    if(nd<0)then; allocate(dmat(0,0)); return; end if
    allocate(dmat(n,nd))
    select case(include)
    case(include_const)
      dmat(:,1)=1.0_dp
    case(include_trend)
      do i=1,n; dmat(i,1)=s+real(i-1,dp); end do
    case(include_both)
      dmat(:,1)=1.0_dp
      do i=1,n; dmat(i,2)=s+real(i-1,dp); end do
    case(include_none)
    end select
  end subroutine build_deterministic

  subroutine lag_embed_univariate(x,m,d,steps,xx,yy,info)
    real(dp),intent(in)::x(:)
    integer,intent(in)::m,d,steps
    real(dp),allocatable,intent(out)::xx(:,:),yy(:)
    integer,intent(out)::info
    integer::nobs,i,j,t
    nobs=size(x)-(m-1)*d-steps
    if(m<1.or.d<1.or.steps<1.or.nobs<1)then
      info=-1; allocate(xx(0,0),yy(0)); return
    end if
    allocate(xx(nobs,m),yy(nobs))
    do i=1,nobs
      t=i+(m-1)*d
      do j=1,m
        xx(i,j)=x(t-(j-1)*d)
      end do
      yy(i)=x(t+steps)
    end do
    info=0
  end subroutine lag_embed_univariate

  subroutine lag_embed_multivariate(y,p,xlag,target,info)
    real(dp),intent(in)::y(:,:)
    integer,intent(in)::p
    real(dp),allocatable,intent(out)::xlag(:,:),target(:,:)
    integer,intent(out)::info
    integer::n,k,nobs,i,j
    n=size(y,1); k=size(y,2); nobs=n-p
    if(p<1.or.nobs<1)then
      info=-1; allocate(xlag(0,0),target(0,0)); return
    end if
    allocate(xlag(nobs,k*p),target(nobs,k)); target=y(p+1:n,:)
    do i=1,nobs
      do j=1,p
        xlag(i,(j-1)*k+1:j*k)=y(p+i-j,:)
      end do
    end do
    info=0
  end subroutine lag_embed_multivariate

  subroutine differences(x,dx)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::dx(:,:)
    allocate(dx(max(0,size(x,1)-1),size(x,2)))
    if(size(x,1)>1)dx=x(2:,:)-x(:size(x,1)-1,:)
  end subroutine differences

  subroutine argsort_real(x,idx)
    real(dp),intent(in)::x(:)
    integer,allocatable,intent(out)::idx(:)
    integer::i,j,key
    allocate(idx(size(x))); idx=[(i,i=1,size(x))]
    do i=2,size(x)
      key=idx(i); j=i-1
      do while(j>=1)
        if(x(idx(j))<x(key))exit
        if(abs(x(idx(j))-x(key))<=epsilon(1.0_dp)*max(1.0_dp,abs(x(key))).and.idx(j)<key)exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=key
    end do
  end subroutine argsort_real

  subroutine sorted_unique(x,u)
    real(dp),intent(in)::x(:)
    real(dp),allocatable,intent(out)::u(:)
    integer,allocatable::idx(:)
    real(dp),allocatable::tmp(:)
    integer::i,n
    if(size(x)==0)then; allocate(u(0)); return; end if
    call argsort_real(x,idx); allocate(tmp(size(x))); n=1; tmp(1)=x(idx(1))
    do i=2,size(x)
      if(abs(x(idx(i))-tmp(n))>epsilon(1.0_dp)*max(1.0_dp,abs(tmp(n))))then; n=n+1; tmp(n)=x(idx(i)); end if
    end do
    allocate(u(n)); u=tmp(1:n)
  end subroutine sorted_unique

  real(dp) function quantile_linear(x,p) result(q)
    real(dp),intent(in)::x(:),p
    integer,allocatable::idx(:)
    real(dp)::h,f
    integer::lo,hi,n
    n=size(x)
    if(n==0)then;q=0.0_dp;return;end if
    call argsort_real(x,idx)
    h=1.0_dp+max(0.0_dp,min(1.0_dp,p))*real(n-1,dp)
    lo=floor(h); hi=ceiling(h); f=h-real(lo,dp)
    q=(1.0_dp-f)*x(idx(lo))+f*x(idx(hi))
  end function quantile_linear

  subroutine percentile_grid(x,trim,ngrid,g)
    real(dp),intent(in)::x(:),trim
    integer,intent(in)::ngrid
    real(dp),allocatable,intent(out)::g(:)
    integer::i
    allocate(g(max(1,ngrid)))
    if(ngrid<=1)then
      g(1)=quantile_linear(x,0.5_dp)
    else
      do i=1,ngrid
        g(i)=quantile_linear(x,trim+(1.0_dp-2.0_dp*trim)*real(i-1,dp)/real(ngrid-1,dp))
      end do
    end if
  end subroutine percentile_grid

  pure elemental real(dp) function sigmoid(x) result(v)
    real(dp),intent(in)::x
    if(x>=0.0_dp)then
      v=1.0_dp/(1.0_dp+exp(-min(x,700.0_dp)))
    else
      v=exp(max(x,-700.0_dp))/(1.0_dp+exp(max(x,-700.0_dp)))
    end if
  end function sigmoid

  pure elemental real(dp) function normal_pdf(x) result(v)
    real(dp),intent(in)::x
    v=exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi_dp)
  end function normal_pdf

  pure elemental real(dp) function normal_cdf(x) result(v)
    real(dp),intent(in)::x
    v=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  subroutine sample_with_replacement(x,out)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::out(:)
    real(dp)::u
    integer::i,j,n
    n=size(x)
    do i=1,size(out)
      call random_number(u); j=min(n,1+int(u*real(n,dp))); out(i)=x(j)
    end do
  end subroutine sample_with_replacement

  pure real(dp) function mean_value(x) result(v)
    real(dp),intent(in)::x(:)
    if(size(x)>0)then;v=sum(x)/real(size(x),dp);else;v=0.0_dp;end if
  end function mean_value

  pure real(dp) function variance_value(x) result(v)
    real(dp),intent(in)::x(:)
    real(dp)::m
    if(size(x)>1)then;m=mean_value(x);v=sum((x-m)**2)/real(size(x)-1,dp);else;v=0.0_dp;end if
  end function variance_value
end module tsdyn_utils
