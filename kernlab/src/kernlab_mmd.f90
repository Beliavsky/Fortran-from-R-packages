! SPDX-License-Identifier: GPL-2.0-only
module kernlab_mmd
  use kernlab_kinds
  use kernlab_types
  use kernlab_kernels, only: kernel_matrix
  use kernlab_linalg, only: quantile_value
  implicit none
  private
  public :: kmmd, kmmd_from_kernels
contains
  subroutine kmmd(x,y,kernel,result,alpha,bootstrap,ntimes,frac)
    real(dp),intent(in)::x(:,:),y(:,:)
    type(kernel_spec),intent(in)::kernel
    type(mmd_result),intent(out)::result
    real(dp),intent(in),optional::alpha,frac
    logical,intent(in),optional::bootstrap
    integer,intent(in),optional::ntimes
    real(dp),allocatable::kxx(:,:),kyy(:,:),kxy(:,:)
    integer::st
    call kernel_matrix(kernel,x,kxx,st)
    if(st/=KL_SUCCESS)then;result%status=st;return;end if
    call kernel_matrix(kernel,y,kyy,st)
    if(st/=KL_SUCCESS)then;result%status=st;return;end if
    call kernel_matrix(kernel,x,kxy,st,y)
    if(st/=KL_SUCCESS)then;result%status=st;return;end if
    call kmmd_from_kernels(kxx,kyy,kxy,result,alpha,bootstrap,ntimes,frac)
  end subroutine kmmd

  subroutine kmmd_from_kernels(kxx,kyy,kxy,result,alpha,bootstrap,ntimes,frac)
    real(dp),intent(in)::kxx(:,:),kyy(:,:),kxy(:,:)
    type(mmd_result),intent(out)::result
    real(dp),intent(in),optional::alpha,frac
    logical,intent(in),optional::bootstrap
    integer,intent(in),optional::ntimes
    integer::m,n,mm,i,j,bn,t,ns
    real(dp)::a,sumxx,sumyy,sumxy,rm,hu,ff
    real(dp),allocatable::boots(:)
    logical::doboot
    result%status=KL_INVALID_ARGUMENT
    m=size(kxx,1);n=size(kyy,1)
    if(size(kxx,2)/=m.or.size(kyy,2)/=n.or.size(kxy,1)/=m.or. &
       size(kxy,2)/=n.or.min(m,n)<2)return
    a=0.05_dp;if(present(alpha))a=alpha;mm=min(m,n)
    sumxx=sum(kxx);sumyy=sum(kyy);sumxy=sum(kxy)
    result%mmd1=sqrt(max(0.0_dp,sumxx/real(m*m,dp)+sumyy/real(n*n,dp) &
      -2.0_dp*sumxy/real(m*n,dp)))
    hu=0.0_dp
    do i=1,mm
      do j=1,mm
        if(i/=j)hu=hu+kxx(i,j)+kyy(i,j)-kxy(i,j)-kxy(j,i)
      end do
    end do
    result%mmd3=hu/real(mm*(mm-1),dp)
    rm=0.0_dp
    do i=1,mm;rm=max(rm,kxx(i,i),kyy(i,i));end do
    result%rademacher_bound=2.0_dp*sqrt(rm/real(mm,dp)) &
      +sqrt(4.0_dp*rm*log(1.0_dp/a)/real(mm,dp))
    result%reject_rademacher=result%mmd1>result%rademacher_bound
    doboot=.false.;if(present(bootstrap))doboot=bootstrap
    if(doboot)then
      bn=100;if(present(ntimes))bn=max(10,ntimes)
      ff=1.0_dp;if(present(frac))ff=max(0.1_dp,min(1.0_dp,frac))
      ns=max(2,ceiling(ff*real(mm,dp)));allocate(boots(bn))
      do t=1,bn
        hu=0.0_dp
        do i=1,ns
          do j=1,ns
            if(i==j)cycle
            hu=hu+bootstrap_term(i,j,t,m,n,kxx,kyy,kxy)
          end do
        end do
        boots(t)=hu/real(ns*(ns-1),dp)
      end do
      result%bootstrap_bound=quantile_value(boots,1.0_dp-a)
      result%reject_bootstrap=result%mmd3>result%bootstrap_bound
    end if
    result%status=KL_SUCCESS
  end subroutine kmmd_from_kernels

  pure real(dp) function bootstrap_term(i,j,t,m,n,kxx,kyy,kxy)
    integer,intent(in)::i,j,t,m,n
    real(dp),intent(in)::kxx(:,:),kyy(:,:),kxy(:,:)
    integer::ii,jj,li,lj,total
    total=m+n
    ii=1+mod(37*i+53*t-1,total);jj=1+mod(37*j+97*t-1,total)
    li=1+mod(19*i+29*t-1,total);lj=1+mod(19*j+71*t-1,total)
    bootstrap_term=combined_kernel(ii,jj,m,kxx,kyy,kxy) &
      -combined_kernel(ii,lj,m,kxx,kyy,kxy) &
      -combined_kernel(li,jj,m,kxx,kyy,kxy) &
      +combined_kernel(li,lj,m,kxx,kyy,kxy)
  end function bootstrap_term

  pure real(dp) function combined_kernel(a,b,m,kxx,kyy,kxy)
    integer,intent(in)::a,b,m
    real(dp),intent(in)::kxx(:,:),kyy(:,:),kxy(:,:)
    if(a<=m.and.b<=m)then
      combined_kernel=kxx(a,b)
    else if(a>m.and.b>m)then
      combined_kernel=kyy(a-m,b-m)
    else if(a<=m)then
      combined_kernel=kxy(a,b-m)
    else
      combined_kernel=kxy(b,a-m)
    end if
  end function combined_kernel
end module kernlab_mmd
