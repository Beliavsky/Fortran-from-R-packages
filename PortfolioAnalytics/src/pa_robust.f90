! SPDX-License-Identifier: GPL-3.0-only
module pa_robust
  use pa_kinds, only: dp
  use pa_statistics, only: sample_moments
  implicit none
  private
  public :: winsorize_returns, covariance_shrinkage, robust_covariance_huber

contains

  subroutine winsorize_returns(returns, lower_fraction, upper_fraction, cleaned)
    real(dp), intent(in) :: returns(:,:), lower_fraction, upper_fraction
    real(dp), intent(out) :: cleaned(:,:)
    real(dp), allocatable :: x(:)
    real(dp) :: lo, hi
    integer :: nobs,nassets,i,j,ilo,ihi
    nobs=size(returns,1)
    nassets=size(returns,2)
    allocate(x(nobs))
    cleaned=returns
    ilo=max(1,min(nobs,1+int(lower_fraction*real(nobs-1,dp))))
    ihi=max(1,min(nobs,1+int(upper_fraction*real(nobs-1,dp))))
    do j=1,nassets
      x=returns(:,j)
      call insertion_sort(x)
      lo=x(ilo)
      hi=x(ihi)
      do i=1,nobs
        cleaned(i,j)=min(max(returns(i,j),lo),hi)
      end do
    end do
  end subroutine winsorize_returns

  subroutine covariance_shrinkage(returns, intensity, mu, sigma)
    real(dp), intent(in) :: returns(:,:), intensity
    real(dp), intent(out) :: mu(:),sigma(:,:)
    real(dp),allocatable::sample(:,:),target(:,:)
    real(dp)::lambda,avgvar
    integer::n,i
    n=size(returns,2)
    allocate(sample(n,n),target(n,n))
    call sample_moments(returns,mu,sample)
    avgvar=sum([(sample(i,i),i=1,n)])/real(n,dp)
    target=0.0_dp
    do i=1,n
      target(i,i)=avgvar
    end do
    lambda=min(max(intensity,0.0_dp),1.0_dp)
    sigma=(1.0_dp-lambda)*sample+lambda*target
  end subroutine covariance_shrinkage

  subroutine robust_covariance_huber(returns, tuning, max_iterations, mu, sigma)
    real(dp), intent(in) :: returns(:,:), tuning
    integer, intent(in) :: max_iterations
    real(dp), intent(out) :: mu(:),sigma(:,:)
    real(dp),allocatable::weights(:),old_mu(:),scale(:),res(:,:),new_mu(:),new_sigma(:,:)
    real(dp)::normv,total
    integer::nobs,nassets,i,j,iter
    nobs=size(returns,1)
    nassets=size(returns,2)
    allocate(weights(nobs),old_mu(nassets),scale(nassets),res(nobs,nassets), &
             new_mu(nassets),new_sigma(nassets,nassets))
    call sample_moments(returns,mu,sigma)
    do iter=1,max_iterations
      old_mu=mu
      do j=1,nassets
        scale(j)=sqrt(max(sigma(j,j),1.0e-16_dp))
      end do
      do i=1,nobs
        res(i,:)=(returns(i,:)-mu)/scale
        normv=sqrt(sum(res(i,:)**2))
        weights(i)=min(1.0_dp,tuning/max(normv,1.0e-16_dp))
      end do
      total=sum(weights)
      new_mu=matmul(transpose(returns),weights)/total
      new_sigma=0.0_dp
      do i=1,nobs
        new_sigma=new_sigma+weights(i)*outer(returns(i,:)-new_mu,returns(i,:)-new_mu)
      end do
      new_sigma=new_sigma/max(total-1.0_dp,1.0_dp)
      mu=new_mu
      sigma=0.5_dp*(new_sigma+transpose(new_sigma))
      if(maxval(abs(mu-old_mu))<1.0e-10_dp) exit
    end do
  end subroutine robust_covariance_huber

  pure function outer(x,y) result(a)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::a(size(x),size(y))
    integer::i
    do i=1,size(x)
      a(i,:)=x(i)*y
    end do
  end function outer

  subroutine insertion_sort(x)
    real(dp),intent(inout)::x(:)
    real(dp)::key
    integer::i,j
    do i=2,size(x)
      key=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=key) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine insertion_sort

end module pa_robust
