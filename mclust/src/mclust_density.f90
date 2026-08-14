! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_density
  use mclust_kinds, only : dp, pi_dp
  use mclust_types, only : mclust_fit
  use mclust_math, only : dmvnorm, logsumexp
  implicit none
  private
  public :: component_log_density, mixture_log_density, cdf_mclust_1d
  public :: quantile_mclust_1d, hdr_levels

contains

  subroutine component_log_density(fit,x,logdens,status)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::logdens(:,:)
    integer,intent(out),optional::status
    real(dp),allocatable::ld(:)
    integer::k,info
    allocate(logdens(size(x,1),fit%g),ld(size(x,1)))
    do k=1,fit%g
      call dmvnorm(x,fit%mean(:,k),fit%sigma(:,:,k),ld,info)
      if(info/=0) then; if(present(status)) status=10+k; return; end if
      logdens(:,k)=ld
    end do
    if(present(status)) status=0
  end subroutine component_log_density

  subroutine mixture_log_density(fit,x,logdens,status)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::logdens(:)
    integer,intent(out),optional::status
    real(dp),allocatable::cd(:,:),row(:)
    integer::i,k,info
    call component_log_density(fit,x,cd,info)
    if(info/=0) then; allocate(logdens(0)); if(present(status)) status=info; return; end if
    allocate(logdens(size(x,1)),row(fit%g))
    do i=1,size(x,1)
      do k=1,fit%g
        row(k)=cd(i,k)+log(max(fit%pro(k),tiny(1.0_dp)))
      end do
      logdens(i)=logsumexp(row)
    end do
    if(present(status)) status=0
  end subroutine mixture_log_density

  pure real(dp) function normal_cdf(x) result(p)
    real(dp),intent(in)::x
    p=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  subroutine cdf_mclust_1d(fit,x,cdf,status)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::cdf(:)
    integer,intent(out),optional::status
    integer::i,k
    real(dp)::sd
    if(fit%d/=1 .or. size(cdf)/=size(x)) then; if(present(status)) status=-1; return; end if
    cdf=0.0_dp
    do k=1,fit%g
      sd=sqrt(fit%sigma(1,1,k))
      do i=1,size(x)
        cdf(i)=cdf(i)+fit%pro(k)*normal_cdf((x(i)-fit%mean(1,k))/sd)
      end do
    end do
    cdf=max(0.0_dp,min(1.0_dp,cdf))
    if(present(status)) status=0
  end subroutine cdf_mclust_1d

  subroutine quantile_mclust_1d(fit,p,q,status,tol,max_iter)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::p(:)
    real(dp),intent(out)::q(:)
    integer,intent(out),optional::status
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::max_iter
    real(dp)::lo,hi,mid,pm,epsq,sdmax
    real(dp)::tmpx(1),tmpcdf(1)
    integer::i,it,nit
    if(fit%d/=1 .or. size(q)/=size(p) .or. any(p<0.0_dp) .or. any(p>1.0_dp)) then
      if(present(status)) status=-1; return
    end if
    epsq=1e-10_dp; if(present(tol)) epsq=tol
    nit=200; if(present(max_iter)) nit=max_iter
    sdmax=maxval(sqrt(fit%sigma(1,1,:)))
    lo=minval(fit%mean(1,:))-12.0_dp*sdmax
    hi=maxval(fit%mean(1,:))+12.0_dp*sdmax
    do i=1,size(p)
      if(p(i)<=0.0_dp) then; q(i)=-huge(1.0_dp); cycle; end if
      if(p(i)>=1.0_dp) then; q(i)=huge(1.0_dp); cycle; end if
      do it=1,nit
        mid=0.5_dp*(lo+hi); tmpx(1)=mid
        call cdf_mclust_1d(fit,tmpx,tmpcdf)
        pm=tmpcdf(1)
        if(pm<p(i)) then; lo=mid; else; hi=mid; end if
        if(abs(hi-lo)<=epsq*(1.0_dp+abs(mid))) exit
      end do
      q(i)=0.5_dp*(lo+hi)
      lo=minval(fit%mean(1,:))-12.0_dp*sdmax
      hi=maxval(fit%mean(1,:))+12.0_dp*sdmax
    end do
    if(present(status)) status=0
  end subroutine quantile_mclust_1d

  subroutine hdr_levels(density,prob,levels)
    real(dp),intent(in)::density(:),prob(:)
    real(dp),intent(out)::levels(:)
    real(dp),allocatable::s(:)
    real(dp)::pos,w
    integer::i,j,n
    n=size(density); allocate(s(n)); s=density; call sort_real(s)
    do i=1,size(prob)
      pos=(1.0_dp-max(0.0_dp,min(1.0_dp,prob(i))))*real(max(0,n-1),dp)+1.0_dp
      j=max(1,min(n,int(floor(pos)))); w=pos-real(j,dp)
      if(j<n) then; levels(i)=(1.0_dp-w)*s(j)+w*s(j+1); else; levels(i)=s(n); end if
    end do
  end subroutine hdr_levels

  subroutine sort_real(a)
    real(dp),intent(inout)::a(:)
    integer::i,j
    real(dp)::v
    do i=2,size(a)
      v=a(i); j=i-1
      do while(j>=1)
        if(a(j)<=v) exit
        a(j+1)=a(j); j=j-1
      end do
      a(j+1)=v
    end do
  end subroutine sort_real
end module mclust_density
