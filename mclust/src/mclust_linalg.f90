! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_linalg
  use mclust_kinds, only : dp
  use r_linalg, only : shared_symmetric_eigen => symmetric_eigen
  implicit none
  private
  public :: symmetric_eigen, inverse_sqrt_symmetric, sample_covariance
contains
  subroutine symmetric_eigen(a,evalues,evec,status,descending)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::evalues(:),evec(:,:)
    integer,intent(out),optional::status
    logical,intent(in),optional::descending
    integer::n,info
    logical::desc
    n=size(a,1); desc=.true.; if(present(descending)) desc=descending
    call shared_symmetric_eigen(a,evalues,evec,info,descending=desc)
    if(present(status)) status=info
  end subroutine symmetric_eigen

  subroutine inverse_sqrt_symmetric(a,invsqrt,inv,status,tol)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::invsqrt(:,:)
    real(dp),allocatable,intent(out),optional::inv(:,:)
    integer,intent(out),optional::status
    real(dp),intent(in),optional::tol
    real(dp),allocatable::ev(:),v(:,:),tmp(:,:)
    real(dp)::thr
    integer::i,info,n
    n=size(a,1); call symmetric_eigen(a,ev,v,info)
    if(info/=0 .or. ev(1)<=0.0_dp) then
      allocate(invsqrt(0,0)); if(present(inv)) allocate(inv(0,0)); if(present(status)) status=merge(info,-1,info/=0); return
    end if
    thr=sqrt(epsilon(1.0_dp)); if(present(tol)) thr=tol
    allocate(tmp(n,n),invsqrt(n,n)); tmp=v
    do i=1,n
      if(ev(i)>max(thr*ev(1),0.0_dp)) then; tmp(:,i)=tmp(:,i)/sqrt(ev(i)); else; tmp(:,i)=0.0_dp; end if
    end do
    invsqrt=matmul(tmp,transpose(v))
    if(present(inv)) inv=matmul(invsqrt,invsqrt)
    if(present(status)) status=0
  end subroutine inverse_sqrt_symmetric

  subroutine sample_covariance(x,mean,cov,population)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::mean(:),cov(:,:)
    logical,intent(in),optional::population
    logical::pop
    integer::i,n
    real(dp)::den
    n=size(x,1); pop=.true.; if(present(population)) pop=population
    mean=sum(x,dim=1)/real(n,dp); cov=0.0_dp
    do i=1,n; cov=cov+outer(x(i,:)-mean,x(i,:)-mean); end do
    den=real(merge(n,max(1,n-1),pop),dp); cov=cov/den
  end subroutine sample_covariance
  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:); real(dp)::c(size(a),size(b)); integer::i
    do i=1,size(a); c(i,:)=a(i)*b; end do
  end function outer
end module mclust_linalg
