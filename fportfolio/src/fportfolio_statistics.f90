! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_statistics
  use fportfolio_kinds, only: dp
  implicit none
  private
  public :: sample_estimator, ewma_estimator, lpm_estimator, slpm_estimator, &
            spearman_estimator, kendall_estimator, shrinkage_estimator, &
            column_means, sample_covariance, quantile_value, skewness_value, kurtosis_value
contains
  pure function column_means(x) result(mu)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: mu(size(x,2))
    mu=sum(x,dim=1)/real(size(x,1),dp)
  end function column_means

  pure function sample_covariance(x) result(sigma)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: sigma(size(x,2),size(x,2))
    real(dp) :: mu(size(x,2)),xc(size(x,1),size(x,2))
    integer :: n
    n=size(x,1); mu=column_means(x)
    xc=x-spread(mu,1,n)
    if (n>1) then
      sigma=matmul(transpose(xc),xc)/real(n-1,dp)
    else
      sigma=0.0_dp
    end if
  end function sample_covariance

  subroutine sample_estimator(x,mu,sigma)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    allocate(mu(size(x,2)),sigma(size(x,2),size(x,2)))
    mu=column_means(x); sigma=sample_covariance(x)
  end subroutine sample_estimator

  subroutine ewma_estimator(x,lambda,mu,sigma)
    real(dp), intent(in) :: x(:,:),lambda
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    real(dp), allocatable :: d(:)
    integer :: i,n,p
    n=size(x,1); p=size(x,2)
    allocate(mu(p),sigma(p,p),d(p)); mu=x(1,:); sigma=0.0_dp
    do i=2,n
      d=x(i,:)-mu
      sigma=lambda*sigma+(1.0_dp-lambda)*spread(d,2,p)*spread(d,1,p)
      mu=lambda*mu+(1.0_dp-lambda)*x(i,:)
    end do
    sigma=0.5_dp*(sigma+transpose(sigma))
  end subroutine ewma_estimator

  subroutine lpm_estimator(x,tau,order,mu,sigma)
    real(dp), intent(in) :: x(:,:),tau(:),order
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    real(dp), allocatable :: z(:,:)
    integer :: n,p
    n=size(x,1);p=size(x,2);allocate(mu(p),sigma(p,p),z(n,p))
    mu=column_means(x)
    z=max(0.0_dp,spread(tau,1,n)-x)**(0.5_dp*order)
    sigma=matmul(transpose(z),z)/real(n,dp)
  end subroutine lpm_estimator

  subroutine slpm_estimator(x,tau,order,mu,sigma)
    real(dp), intent(in) :: x(:,:),tau(:),order
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    real(dp), allocatable :: z(:,:)
    integer :: n,p
    n=size(x,1);p=size(x,2);allocate(mu(p),sigma(p,p),z(n,p))
    mu=column_means(x)
    z=sign(abs(x-spread(tau,1,n))**(0.5_dp*order),x-spread(tau,1,n))
    sigma=matmul(transpose(z),z)/real(n,dp)
  end subroutine slpm_estimator

  subroutine spearman_estimator(x,mu,sigma)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    real(dp), allocatable :: ranks(:,:),sds(:)
    integer :: j,p,n
    n=size(x,1);p=size(x,2);allocate(mu(p),sigma(p,p),ranks(n,p),sds(p))
    mu=column_means(x)
    do j=1,p
      call average_ranks(x(:,j),ranks(:,j))
      sds(j)=sqrt(max(sum((x(:,j)-mu(j))**2)/real(max(1,n-1),dp),0.0_dp))
    end do
    sigma=sample_covariance(ranks)
    do j=1,p
      sigma(j,j)=1.0_dp
    end do
    call covariance_from_correlation(sigma,sds)
  end subroutine spearman_estimator

  subroutine kendall_estimator(x,mu,sigma)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    real(dp), allocatable :: sds(:)
    integer :: i,j,p,n,a,b
    real(dp) :: concord,tau
    n=size(x,1);p=size(x,2);allocate(mu(p),sigma(p,p),sds(p));mu=column_means(x)
    do j=1,p
      sds(j)=sqrt(max(sum((x(:,j)-mu(j))**2)/real(max(1,n-1),dp),0.0_dp))
    end do
    sigma=0.0_dp
    do i=1,p
      sigma(i,i)=1.0_dp
      do j=i+1,p
        concord=0.0_dp
        do a=1,n-1; do b=a+1,n
          concord=concord+sign(1.0_dp,(x(b,i)-x(a,i))*(x(b,j)-x(a,j)))
          if (abs((x(b,i)-x(a,i))*(x(b,j)-x(a,j))) <= tiny(1.0_dp)) concord=concord-1.0_dp
        end do; end do
        tau=2.0_dp*concord/real(n*(n-1),dp)
        sigma(i,j)=tau;sigma(j,i)=tau
      end do
    end do
    call covariance_from_correlation(sigma,sds)
  end subroutine kendall_estimator

  subroutine shrinkage_estimator(x,mu,sigma,intensity)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mu(:),sigma(:,:)
    real(dp), intent(in), optional :: intensity
    real(dp), allocatable :: s(:,:),target(:,:),xc(:,:)
    real(dp) :: lambda,phi,rho
    integer :: n,p,i
    call sample_estimator(x,mu,s)
    n=size(x,1);p=size(x,2);allocate(target(p,p),xc(n,p),sigma(p,p))
    target=0.0_dp
    do i=1,p; target(i,i)=s(i,i); end do
    if (present(intensity)) then
      lambda=max(0.0_dp,min(1.0_dp,intensity))
    else
      xc=x-spread(mu,1,n);phi=0.0_dp
      do i=1,n
        phi=phi+sum((spread(xc(i,:),2,p)*spread(xc(i,:),1,p)-s)**2)
      end do
      phi=phi/real(n,dp);rho=sum((s-target)**2)
      if (rho>0.0_dp) then;lambda=max(0.0_dp,min(1.0_dp,phi/(real(n,dp)*rho)));else;lambda=1.0_dp;end if
    end if
    sigma=(1.0_dp-lambda)*s+lambda*target
  end subroutine shrinkage_estimator

  subroutine covariance_from_correlation(cor,sds)
    real(dp), intent(inout) :: cor(:,:)
    real(dp), intent(in) :: sds(:)
    integer :: i,j
    do i=1,size(cor,1); do j=1,size(cor,2); cor(i,j)=cor(i,j)*sds(i)*sds(j); end do; end do
  end subroutine covariance_from_correlation

  subroutine average_ranks(x,r)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: r(:)
    integer :: idx(size(x)),i,j,k,n,tmp
    n=size(x);idx=[(i,i=1,n)]
    do i=2,n
      tmp=idx(i)
      j=i-1
      do while(j>=1)
        if(x(idx(j))<=x(tmp))exit
        idx(j+1)=idx(j)
        j=j-1
      end do
      idx(j+1)=tmp
    end do
    i=1
    do while(i<=n)
      j=i
      do while(j<n)
        if(abs(x(idx(j+1))-x(idx(i)))>tiny(1.0_dp))exit
        j=j+1
      end do
      do k=i,j
        r(idx(k))=0.5_dp*real(i+j,dp)
      end do
      i=j+1
    end do
  end subroutine average_ranks

  real(dp) function quantile_value(x,p) result(q)
    real(dp), intent(in) :: x(:),p
    real(dp), allocatable :: y(:)
    real(dp) :: h
    integer :: i,j,n
    y=x;n=size(y)
    call sort_values(y)
    if(n==0)then;q=0.0_dp;return;end if
    h=1.0_dp+(real(n-1,dp))*max(0.0_dp,min(1.0_dp,p));i=floor(h);j=min(n,i+1)
    q=y(i)+(h-real(i,dp))*(y(j)-y(i))
  end function quantile_value

  subroutine sort_values(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_values

  pure real(dp) function skewness_value(x) result(v)
    real(dp), intent(in)::x(:)
    real(dp)::m,s
    m=sum(x)/real(size(x),dp);s=sqrt(sum((x-m)**2)/real(max(1,size(x)-1),dp))
    if(s>0.0_dp)then;v=sum(((x-m)/s)**3)/real(size(x),dp);else;v=0.0_dp;end if
  end function skewness_value

  pure real(dp) function kurtosis_value(x) result(v)
    real(dp), intent(in)::x(:)
    real(dp)::m,s
    m=sum(x)/real(size(x),dp);s=sqrt(sum((x-m)**2)/real(max(1,size(x)-1),dp))
    if(s>0.0_dp)then;v=sum(((x-m)/s)**4)/real(size(x),dp);else;v=3.0_dp;end if
  end function kurtosis_value
end module fportfolio_statistics
