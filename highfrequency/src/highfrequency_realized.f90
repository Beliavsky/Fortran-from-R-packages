! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_realized
  use highfrequency_kinds, only: dp, pi
  use highfrequency_stats, only: median3, mu_abs_normal
  use highfrequency_linalg, only: covariance_to_correlation, make_psd, symmetrize
  use highfrequency_data, only: make_returns
  implicit none
  private

  public :: realized_variance, realized_covariance, realized_skewness
  public :: realized_kurtosis, realized_semivariance, realized_semicovariance
  public :: bipower_variation, bipower_covariance, realized_quarticity
  public :: tripower_quarticity, quadpower_variation, minimum_realized_variance
  public :: median_realized_variance, minimum_realized_quarticity
  public :: median_realized_quarticity, multipower_variation
  public :: threshold_covariance, realized_kernel_covariance
  public :: two_scale_variance, two_scale_covariance
  public :: hayashi_yoshida_covariance, realized_beta, noise_variance
  public :: modulated_realized_covariance
  public :: preaveraged_covariance, average_realized_covariance

  interface realized_variance
    module procedure realized_variance_vector
    module procedure realized_variance_matrix
  end interface realized_variance

contains

  pure real(dp) function realized_variance_vector(r) result(value)
    real(dp), intent(in) :: r(:)
    value = sum(r*r)
  end function realized_variance_vector

  pure function realized_variance_matrix(r) result(value)
    real(dp), intent(in) :: r(:,:)
    real(dp) :: value(size(r,2))
    integer :: j
    do j = 1, size(r,2)
      value(j) = sum(r(:,j)*r(:,j))
    end do
  end function realized_variance_matrix

  pure function realized_covariance(r, correlation) result(value)
    real(dp), intent(in) :: r(:,:)
    logical, intent(in), optional :: correlation
    real(dp) :: value(size(r,2),size(r,2))
    real(dp) :: cov(size(r,2),size(r,2))
    cov = matmul(transpose(r),r)
    call symmetrize(cov)
    if (present(correlation)) then
      if (correlation) then
        call covariance_to_correlation(cov,value)
      else
        value=cov
      end if
    else
      value=cov
    end if
  end function realized_covariance

  pure real(dp) function realized_skewness(r) result(value)
    real(dp), intent(in) :: r(:)
    real(dp) :: rv
    rv=sum(r*r)
    if (rv <= 0.0_dp) then
      value=0.0_dp
    else
      value=sqrt(real(size(r),dp))*sum(r**3)/(rv**1.5_dp)
    end if
  end function realized_skewness

  pure real(dp) function realized_kurtosis(r) result(value)
    real(dp), intent(in) :: r(:)
    real(dp) :: rv
    rv=sum(r*r)
    if (rv <= 0.0_dp) then
      value=0.0_dp
    else
      value=real(size(r),dp)*sum(r**4)/(rv*rv)
    end if
  end function realized_kurtosis

  pure function realized_semivariance(r) result(value)
    real(dp), intent(in) :: r(:)
    real(dp) :: value(2)
    value=0.0_dp
    value(1)=sum(merge(r*r,0.0_dp,r<0.0_dp))
    value(2)=sum(merge(r*r,0.0_dp,r>0.0_dp))
  end function realized_semivariance

  pure function realized_semicovariance(r, side) result(value)
    real(dp), intent(in) :: r(:,:)
    integer, intent(in), optional :: side
    real(dp) :: value(size(r,2),size(r,2))
    real(dp), allocatable :: rr(:,:)
    integer :: s, i, j
    s=-1
    if(present(side)) s=side
    allocate(rr(size(r,1),size(r,2)))
    rr=0.0_dp
    if(s<0)then
      where(r<0.0_dp) rr=r
    else if(s>0)then
      where(r>0.0_dp) rr=r
    else
      rr=r
    end if
    value=matmul(transpose(rr),rr)
    do i=1,size(value,1)
      do j=i+1,size(value,2)
        value(i,j)=0.5_dp*(value(i,j)+value(j,i))
        value(j,i)=value(i,j)
      end do
    end do
  end function realized_semicovariance

  pure real(dp) function bipower_variation(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: n
    n=size(r)
    if(n<2)then
      value=0.0_dp
    else
      value=(pi/2.0_dp)*real(n,dp)/real(n-1,dp)*sum(abs(r(1:n-1))*abs(r(2:n)))
    end if
  end function bipower_variation

  pure function bipower_covariance(r, correlation) result(value)
    real(dp), intent(in) :: r(:,:)
    logical, intent(in), optional :: correlation
    real(dp) :: value(size(r,2),size(r,2))
    real(dp) :: cov(size(r,2),size(r,2))
    real(dp), allocatable :: plus(:), minus(:)
    integer :: i,j,n
    n=size(r,1)
    cov=0.0_dp
    do i=1,size(r,2)
      cov(i,i)=bipower_variation(r(:,i))
      do j=1,i-1
        allocate(plus(n),minus(n))
        plus=r(:,i)+r(:,j)
        minus=r(:,i)-r(:,j)
        cov(i,j)=0.25_dp*(bipower_variation(plus)-bipower_variation(minus))
        cov(j,i)=cov(i,j)
        deallocate(plus,minus)
      end do
    end do
    if(present(correlation))then
      if(correlation)then
        call covariance_to_correlation(cov,value)
      else
        value=cov
      end if
    else
      value=cov
    end if
  end function bipower_covariance

  pure real(dp) function realized_quarticity(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: n
    n=size(r)
    value=real(n,dp)/3.0_dp*sum(r**4)
  end function realized_quarticity

  pure real(dp) function tripower_quarticity(r) result(value)
    real(dp), intent(in) :: r(:)
    real(dp) :: c
    integer :: i,n
    n=size(r)
    if(n<3)then
      value=0.0_dp
      return
    end if
    c=(gamma(0.5_dp)/(2.0_dp**(2.0_dp/3.0_dp)*gamma(7.0_dp/6.0_dp)))**3
    value=0.0_dp
    do i=1,n-2
      value=value+(abs(r(i)*r(i+1)*r(i+2)))**(4.0_dp/3.0_dp)
    end do
    value=real(n,dp)*real(n,dp)/real(n-2,dp)*c*value
  end function tripower_quarticity

  pure real(dp) function quadpower_variation(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: i,n
    n=size(r)
    if(n<4)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=1,n-3
      value=value+abs(r(i)*r(i+1)*r(i+2)*r(i+3))
    end do
    value=real(n,dp)*real(n,dp)/real(n-3,dp)*(pi*pi/4.0_dp)*value
  end function quadpower_variation

  pure real(dp) function minimum_realized_variance(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: i,n
    n=size(r)
    if(n<2)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=1,n-1
      value=value+min(abs(r(i)),abs(r(i+1)))**2
    end do
    value=pi/(pi-2.0_dp)*real(n,dp)/real(n-1,dp)*value
  end function minimum_realized_variance

  pure real(dp) function median_realized_variance(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: i,n
    n=size(r)
    if(n<3)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=2,n-1
      value=value+median3(abs(r(i-1)),abs(r(i)),abs(r(i+1)))**2
    end do
    value=pi/(6.0_dp-4.0_dp*sqrt(3.0_dp)+pi)*real(n,dp)/real(n-2,dp)*value
  end function median_realized_variance

  pure real(dp) function minimum_realized_quarticity(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: i,n
    n=size(r)
    if(n<2)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=1,n-1
      value=value+min(abs(r(i)),abs(r(i+1)))**4
    end do
    value=pi*real(n,dp)/(3.0_dp*pi-8.0_dp)*real(n,dp)/real(n-1,dp)*value
  end function minimum_realized_quarticity

  pure real(dp) function median_realized_quarticity(r) result(value)
    real(dp), intent(in) :: r(:)
    integer :: i,n
    n=size(r)
    if(n<3)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=2,n-1
      value=value+median3(abs(r(i-1)),abs(r(i)),abs(r(i+1)))**4
    end do
    value=3.0_dp*pi*real(n,dp)/(9.0_dp*pi+72.0_dp-52.0_dp*sqrt(3.0_dp))* &
      real(n,dp)/real(n-2,dp)*value
  end function median_realized_quarticity

  pure real(dp) function multipower_variation(r,m,p) result(value)
    real(dp), intent(in) :: r(:)
    integer, intent(in) :: m
    real(dp), intent(in) :: p
    real(dp) :: product_value, exponent, normalization
    integer :: i,j,n
    n=size(r)
    if(m<1 .or. n<m .or. p<=0.0_dp)then
      value=0.0_dp
      return
    end if
    exponent=p/real(m,dp)
    normalization=mu_abs_normal(exponent)**(-real(m,dp))
    value=0.0_dp
    do i=1,n-m+1
      product_value=1.0_dp
      do j=0,m-1
        product_value=product_value*abs(r(i+j))**exponent
      end do
      value=value+product_value
    end do
    value=real(n,dp)**(0.5_dp*p-1.0_dp)*real(n,dp)/real(n-m+1,dp)*normalization*value
  end function multipower_variation

  function threshold_covariance(r, correlation) result(value)
    real(dp), intent(in) :: r(:,:)
    logical, intent(in), optional :: correlation
    real(dp) :: value(size(r,2),size(r,2))
    real(dp), allocatable :: truncated(:,:)
    real(dp) :: threshold(size(r,2)), bp
    integer :: i,j,n
    n=size(r,1)
    allocate(truncated(n,size(r,2)))
    truncated=r
    do j=1,size(r,2)
      bp=bipower_variation(r(:,j))
      threshold(j)=3.0_dp*sqrt(max(0.0_dp,bp))*(1.0_dp/real(max(1,n),dp))**0.49_dp
      do i=1,n
        if(abs(truncated(i,j))>threshold(j)) truncated(i,j)=0.0_dp
      end do
    end do
    value=realized_covariance(truncated,correlation)
  end function threshold_covariance

  pure real(dp) function kernel_weight(kind,x) result(value)
    character(len=*), intent(in) :: kind
    real(dp), intent(in) :: x
    real(dp) :: z
    z=abs(x)
    select case(trim(adjustl(kind)))
    case('bartlett','Bartlett')
      value=max(0.0_dp,1.0_dp-z)
    case('parzen','Parzen')
      if(z<=0.5_dp)then
        value=1.0_dp-6.0_dp*z*z+6.0_dp*z**3
      else if(z<=1.0_dp)then
        value=2.0_dp*(1.0_dp-z)**3
      else
        value=0.0_dp
      end if
    case('tukey-hanning','TukeyHanning')
      if(z<=1.0_dp)then
        value=0.5_dp*(1.0_dp+cos(pi*z))
      else
        value=0.0_dp
      end if
    case('flat-top','FlatTop')
      if(z<=0.5_dp)then
        value=1.0_dp
      else if(z<=1.0_dp)then
        value=2.0_dp*(1.0_dp-z)
      else
        value=0.0_dp
      end if
    case default
      value=max(0.0_dp,1.0_dp-z)
    end select
  end function kernel_weight

  function realized_kernel_covariance(r, bandwidth, kernel, correlation, force_psd) result(value)
    real(dp), intent(in) :: r(:,:)
    integer, intent(in) :: bandwidth
    character(len=*), intent(in), optional :: kernel
    logical, intent(in), optional :: correlation, force_psd
    real(dp) :: value(size(r,2),size(r,2))
    real(dp) :: cov(size(r,2),size(r,2)), gamma_h(size(r,2),size(r,2))
    real(dp) :: weight
    character(len=32) :: kind
    logical :: cor_flag, psd_flag, ok
    integer :: h,n,i,a,b
    n=size(r,1)
    kind='parzen'
    if(present(kernel)) kind=kernel
    cor_flag=.false.
    if(present(correlation)) cor_flag=correlation
    psd_flag=.false.
    if(present(force_psd)) psd_flag=force_psd
    cov=matmul(transpose(r),r)
    do h=1,min(max(0,bandwidth),n-1)
      gamma_h=0.0_dp
      do i=1,n-h
        do a=1,size(r,2)
          do b=1,size(r,2)
            gamma_h(a,b)=gamma_h(a,b)+r(i+h,a)*r(i,b)
          end do
        end do
      end do
      weight=kernel_weight(kind,real(h,dp)/real(bandwidth+1,dp))
      cov=cov+weight*(gamma_h+transpose(gamma_h))
    end do
    call symmetrize(cov)
    if(psd_flag)then
      call make_psd(cov,value,0.0_dp,ok)
      cov=value
    end if
    if(cor_flag)then
      call covariance_to_correlation(cov,value)
    else
      value=cov
    end if
  end function realized_kernel_covariance

  pure real(dp) function subsampled_variance(log_price,step) result(value)
    real(dp), intent(in) :: log_price(:)
    integer,intent(in)::step
    integer::i,n
    n=size(log_price)
    if(step<1 .or. n<=step)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=1,n-step
      value=value+(log_price(i+step)-log_price(i))**2
    end do
    value=value/real(step,dp)
  end function subsampled_variance

  pure real(dp) function two_scale_variance(prices,k,j) result(value)
    real(dp),intent(in)::prices(:)
    integer,intent(in)::k,j
    real(dp),allocatable::lp(:)
    real(dp)::nk,nj,slow,fast,denom
    integer::n
    n=size(prices)
    if(n<=max(k,j) .or. min(k,j)<1 .or. any(prices<=0.0_dp))then
      value=0.0_dp
      return
    end if
    allocate(lp(n)); lp=log(prices)
    nk=real(n-k,dp)/real(k,dp)
    nj=real(n-j,dp)/real(j,dp)
    slow=subsampled_variance(lp,k)
    fast=subsampled_variance(lp,j)
    denom=1.0_dp-nk/nj
    if(abs(denom)<=epsilon(1.0_dp))then
      value=0.0_dp
    else
      value=(slow-(nk/nj)*fast)/denom
    end if
  end function two_scale_variance

  pure real(dp) function lagged_cross(logx,logy,step) result(value)
    real(dp),intent(in)::logx(:),logy(:)
    integer,intent(in)::step
    integer::i,n
    n=min(size(logx),size(logy))
    if(step<1 .or. n<=step)then
      value=0.0_dp
      return
    end if
    value=0.0_dp
    do i=1,n-step
      value=value+(logx(i+step)-logx(i))*(logy(i+step)-logy(i))
    end do
    value=value/real(step,dp)
  end function lagged_cross

  pure real(dp) function two_scale_covariance(px,py,k,j) result(value)
    real(dp),intent(in)::px(:),py(:)
    integer,intent(in)::k,j
    real(dp),allocatable::lx(:),ly(:)
    real(dp)::nk,nj,slow,fast,denom
    integer::n
    n=min(size(px),size(py))
    if(n<=max(k,j) .or. min(k,j)<1 .or. any(px(:n)<=0.0_dp) .or. any(py(:n)<=0.0_dp))then
      value=0.0_dp
      return
    end if
    allocate(lx(n),ly(n)); lx=log(px(:n));ly=log(py(:n))
    nk=real(n-k,dp)/real(k,dp)
    nj=real(n-j,dp)/real(j,dp)
    slow=lagged_cross(lx,ly,k)
    fast=lagged_cross(lx,ly,j)
    denom=1.0_dp-nk/nj
    if(abs(denom)<=epsilon(1.0_dp))then
      value=0.0_dp
    else
      value=(slow-(nk/nj)*fast)/denom
    end if
  end function two_scale_covariance

  pure real(dp) function hayashi_yoshida_covariance(t1,p1,t2,p2) result(value)
    integer,intent(in)::t1(:),t2(:)
    real(dp),intent(in)::p1(:),p2(:)
    real(dp),allocatable::r1(:),r2(:)
    integer::i,j
    value=0.0_dp
    if(size(p1)<2 .or. size(p2)<2) return
    if(any(p1<=0.0_dp) .or. any(p2<=0.0_dp)) return
    r1=make_returns(p1)
    r2=make_returns(p2)
    do i=1,size(r1)
      do j=1,size(r2)
        if(max(t1(i),t2(j))<min(t1(i+1),t2(j+1))) value=value+r1(i)*r2(j)
      end do
    end do
  end function hayashi_yoshida_covariance

  pure real(dp) function realized_beta(asset,index_returns) result(value)
    real(dp),intent(in)::asset(:),index_returns(:)
    integer::n
    real(dp)::denom
    n=min(size(asset),size(index_returns))
    denom=sum(index_returns(:n)**2)
    if(denom<=0.0_dp)then
      value=0.0_dp
    else
      value=sum(asset(:n)*index_returns(:n))/denom
    end if
  end function realized_beta

  pure real(dp) function noise_variance(returns) result(value)
    real(dp),intent(in)::returns(:)
    integer::n
    n=size(returns)
    if(n<2)then
      value=0.0_dp
    else
      value=-sum(returns(2:n)*returns(1:n-1))/real(n-1,dp)
      value=max(0.0_dp,value)
    end if
  end function noise_variance

  function modulated_realized_covariance(r,theta,correlation,force_psd) result(value)
    real(dp),intent(in)::r(:,:)
    real(dp),intent(in),optional::theta
    logical,intent(in),optional::correlation,force_psd
    real(dp)::value(size(r,2),size(r,2))
    real(dp)::tmp(size(r,2),size(r,2))
    real(dp),allocatable::pre(:,:)
    real(dp)::th
    integer::i,n
    logical::cor_flag,psd_flag,ok
    n=size(r,1)
    th=0.8_dp
    if(present(theta)) th=theta
    allocate(pre(max(1,n-1),size(r,2)))
    if(n<2)then
      pre=0.0_dp
    else
      do i=1,n-1
        pre(i,:)=th*r(i,:)+(1.0_dp-th)*r(i+1,:)
      end do
    end if
    value=matmul(transpose(pre),pre)
    if(n>1) value=value*real(n,dp)/real(n-1,dp)
    cor_flag=.false.;if(present(correlation))cor_flag=correlation
    psd_flag=.false.;if(present(force_psd))psd_flag=force_psd
    if(psd_flag)then
      call make_psd(value,tmp,0.0_dp,ok)
      value=tmp
    end if
    if(cor_flag) then
      call covariance_to_correlation(value,tmp)
      value=tmp
    end if
  end function modulated_realized_covariance


  function preaveraged_covariance(r,window,correlation,force_psd) result(value)
    real(dp),intent(in)::r(:,:)
    integer,intent(in),optional::window
    logical,intent(in),optional::correlation,force_psd
    real(dp)::value(size(r,2),size(r,2))
    real(dp),allocatable::pre(:,:),weights(:),tmp(:,:)
    real(dp)::psi1,psi2
    integer::n,k,i,j
    logical::cor_flag,psd_flag,ok
    n=size(r,1)
    k=max(2,int(sqrt(real(max(1,n),dp))))
    if(present(window))k=max(2,window)
    if(n<k+1)then
      value=0.0_dp
      return
    end if
    allocate(weights(k-1),pre(n-k+1,size(r,2)),tmp(size(r,2),size(r,2)))
    do j=1,k-1
      weights(j)=min(real(j,dp)/real(k,dp),1.0_dp-real(j,dp)/real(k,dp))
    end do
    pre=0.0_dp
    do i=1,n-k+1
      do j=1,k-1
        pre(i,:)=pre(i,:)+weights(j)*r(i+j-1,:)
      end do
    end do
    psi2=sum(weights*weights)/real(k,dp)
    psi1=(weights(1)**2+weights(size(weights))**2)
    if(size(weights)>1)psi1=psi1+sum((weights(2:)-weights(:size(weights)-1))**2)
    psi1=psi1*real(k,dp)
    value=matmul(transpose(pre),pre)/(psi2*real(k,dp))
    value=value-psi1/(2.0_dp*psi2*real(k*k,dp))*matmul(transpose(r),r)
    call symmetrize(value)
    psd_flag=.false.;if(present(force_psd))psd_flag=force_psd
    if(psd_flag)then
      call make_psd(value,tmp,0.0_dp,ok)
      value=tmp
    end if
    cor_flag=.false.;if(present(correlation))cor_flag=correlation
    if(cor_flag)then
      call covariance_to_correlation(value,tmp)
      value=tmp
    end if
  end function preaveraged_covariance

  function average_realized_covariance(r,subsamples,correlation) result(value)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::subsamples
    logical,intent(in),optional::correlation
    real(dp)::value(size(r,2),size(r,2)),tmp(size(r,2),size(r,2))
    real(dp),allocatable::block(:,:)
    integer::n,k,offset,i,j,count,nb
    logical::cor_flag
    n=size(r,1)
    k=max(1,subsamples)
    value=0.0_dp
    count=0
    do offset=1,k
      nb=(n-offset+1)/k
      if(nb<1)cycle
      allocate(block(nb,size(r,2)))
      do i=1,nb
        block(i,:)=0.0_dp
        do j=offset+(i-1)*k,min(n,offset+i*k-1)
          block(i,:)=block(i,:)+r(j,:)
        end do
      end do
      value=value+matmul(transpose(block),block)
      count=count+1
      deallocate(block)
    end do
    if(count>0)value=value/real(count,dp)
    call symmetrize(value)
    cor_flag=.false.;if(present(correlation))cor_flag=correlation
    if(cor_flag)then
      call covariance_to_correlation(value,tmp)
      value=tmp
    end if
  end function average_realized_covariance

end module highfrequency_realized
