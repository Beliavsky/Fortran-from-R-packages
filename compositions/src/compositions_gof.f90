! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_gof
  use compositions_kinds, only: dp, pi
  use bayesm_rng, only: rand_uniform
  implicit none
  private
  public :: kernel_similarity_statistic, kernel_similarity_permutations
  public :: poisson_ks_statistic, sorted_uniforms, poisson_ks_sample
contains
  real(dp) function gaussian_kernel(x,y,bw) result(k)
    real(dp), intent(in) :: x(:),y(:),bw
    real(dp) :: normq,sq
    if(size(x)/=size(y).or.bw<=0.0_dp) error stop 'gaussian_kernel: invalid arguments'
    normq=sum((x-y)**2); sq=2.0_dp*bw*bw*pi
    k=exp(-normq/(2.0_dp*bw*bw))/sqrt(sq**real(size(x),dp))
  end function gaussian_kernel

  real(dp) function kernel_similarity_statistic(x,y,bw) result(stat)
    real(dp), intent(in) :: x(:,:),y(:,:),bw
    real(dp) :: p1,p2,p3
    integer :: i,j,n,m
    if(size(x,2)/=size(y,2)) error stop 'kernel_similarity_statistic: dimension mismatch'
    n=size(x,1); m=size(y,1); p1=0.0_dp; p2=0.0_dp; p3=0.0_dp
    ! This intentionally follows gsiDensityCheck's lower-triangle accumulation.
    do i=1,n
      do j=1,i
        p1=p1+gaussian_kernel(x(i,:),x(j,:),bw)
      end do
    end do
    do i=1,m
      do j=1,n
        p2=p2+gaussian_kernel(y(i,:),x(j,:),bw)
      end do
      do j=1,i
        p3=p3+gaussian_kernel(y(i,:),y(j,:),bw)
      end do
    end do
    p1=p1/real(n*n,dp); p2=p2/real(n*m,dp); p3=p3/real(m*m,dp)
    stat=p2/sqrt(p1*p3)
  end function kernel_similarity_statistic

  function kernel_similarity_permutations(x,y,bw,nreps) result(reps)
    real(dp), intent(in) :: x(:,:),y(:,:),bw
    integer, intent(in) :: nreps
    real(dp) :: reps(nreps)
    real(dp), allocatable :: allx(:,:),a(:,:),b(:,:)
    integer, allocatable :: idx(:)
    integer :: n,m,l,r,i,j,tmp
    n=size(x,1); m=size(y,1); l=n+m
    if(size(x,2)/=size(y,2)) error stop 'kernel_similarity_permutations: dimension mismatch'
    allocate(allx(l,size(x,2)),idx(l),a(n,size(x,2)),b(m,size(x,2)))
    allx(1:n,:)=x; allx(n+1:l,:)=y
    do r=1,nreps
      idx=[(i,i=1,l)]
      do i=1,n
        j=i+int(rand_uniform()*real(l-i+1,dp)); if(j>l) j=l
        tmp=idx(i); idx(i)=idx(j); idx(j)=tmp
      end do
      a=allx(idx(1:n),:); b=allx(idx(n+1:l),:)
      reps(r)=kernel_similarity_statistic(a,b,bw)
    end do
  end function kernel_similarity_permutations

  real(dp) function poisson_ks_statistic(data,ps) result(stat)
    integer, intent(in) :: data(:)
    real(dp), intent(in) :: ps(:)
    integer, allocatable :: counts(:)
    real(dp) :: delta,maxd
    integer :: i,n,v
    n=size(data); allocate(counts(size(ps))); counts=0
    do i=1,n
      v=data(i)
      if(v>=0.and.v<size(ps)) counts(v+1)=counts(v+1)+1
    end do
    delta=0.0_dp; maxd=0.0_dp
    do i=1,size(ps)
      delta=delta+real(n,dp)*ps(i)-real(counts(i),dp)
      maxd=max(maxd,abs(delta))
    end do
    stat=maxd/real(n,dp)
  end function poisson_ks_statistic

  function sorted_uniforms(n) result(data)
    integer, intent(in) :: n
    real(dp) :: data(n),tmp,u
    integer :: i
    tmp=0.0_dp
    do i=1,n
      u=max(rand_uniform(),tiny(1.0_dp)); tmp=tmp-log(u); data(i)=tmp
    end do
    u=max(rand_uniform(),tiny(1.0_dp)); tmp=tmp-log(u)
    data=data/tmp
  end function sorted_uniforms

  function poisson_ks_sample(n,ps,nsamples) result(statistics)
    integer, intent(in) :: n,nsamples
    real(dp), intent(in) :: ps(:)
    real(dp) :: statistics(nsamples),cdf(size(ps)),u(n),tmp
    integer :: r,i,j
    cdf=ps
    do i=2,size(cdf); cdf(i)=cdf(i)+cdf(i-1); end do
    do r=1,nsamples
      u=sorted_uniforms(n); j=1; statistics(r)=0.0_dp
      do i=1,size(cdf)
        do while(j<=n)
          if(u(j)>cdf(i)) exit
          j=j+1
        end do
        tmp=abs(real(j-1,dp)/real(n,dp)-cdf(i))
        statistics(r)=max(statistics(r),tmp)
      end do
    end do
  end function poisson_ks_sample
end module compositions_gof
