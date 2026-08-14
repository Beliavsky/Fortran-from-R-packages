! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_dr
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit
  use mclust_linalg, only : inverse_sqrt_symmetric, symmetric_eigen, sample_covariance
  implicit none
  private
  type, public :: mclust_dr_fit
    integer :: d=0
    integer :: numdir=0
    real(dp) :: lambda=1.0_dp
    real(dp),allocatable :: evalues(:)
    real(dp),allocatable :: basis(:,:)
    real(dp),allocatable :: raw_evectors(:,:)
    real(dp),allocatable :: direction(:,:)
    real(dp),allocatable :: kernel(:,:), kernel_mean(:,:), kernel_cov(:,:)
  end type
  public :: fit_mclust_dr, project_mclust_dr
contains
  subroutine fit_mclust_dr(fit,x,dr,lambda,normalized,status)
    type(mclust_fit),intent(in)::fit
    real(dp),intent(in)::x(:,:)
    type(mclust_dr_fit),intent(out)::dr
    real(dp),intent(in),optional::lambda
    logical,intent(in),optional::normalized
    integer,intent(out),optional::status
    real(dp),allocatable::mu(:),sigma0(:,:),invs(:,:),inv(:,:),s(:,:),mi(:,:),mii(:,:),m(:,:),ev(:),u(:,:),a(:,:),delta(:)
    real(dp)::lam,nrm
    logical::norm
    integer::d,g,j,k,info,nd
    d=size(x,2); g=fit%g; lam=1.0_dp; if(present(lambda)) lam=max(0.0_dp,min(1.0_dp,lambda))
    norm=.true.; if(present(normalized)) norm=normalized
    allocate(mu(d),sigma0(d,d)); call sample_covariance(x,mu,sigma0,population=.true.)
    call inverse_sqrt_symmetric(sigma0,invs,inv,info)
    if(info/=0) then; if(present(status)) status=info; return; end if
    allocate(s(d,d),mi(d,d),mii(d,d),m(d,d),delta(d)); s=0.0_dp; mi=0.0_dp; mii=0.0_dp
    do k=1,g; s=s+fit%pro(k)*fit%sigma(:,:,k); end do
    do k=1,g
      delta=fit%mean(:,k)-mu; mi=mi+fit%pro(k)*outer(delta,delta)
      if(lam<1.0_dp) then
        allocate(a(d,d)); a=matmul(invs,fit%sigma(:,:,k)-s)
        mii=mii+fit%pro(k)*matmul(transpose(a),a); deallocate(a)
      end if
    end do
    a=matmul(invs,mi)
    m=2.0_dp*lam*matmul(transpose(a),a)+2.0_dp*(1.0_dp-lam)*mii
    ! generalized eigensystem of M with respect to Sigma: invsqrt*M*invsqrt
    a=matmul(transpose(invs),matmul(m,invs))
    call symmetric_eigen(a,ev,u,info)
    if(info/=0) then; if(present(status)) status=info; return; end if
    u=matmul(invs,u)
    ev=max(ev,0.0_dp); nd=min(d,count(ev>sqrt(epsilon(1.0_dp))))
    if(nd<1) nd=1
    if(norm) then
      do j=1,nd; nrm=sqrt(dot_product(u(:,j),u(:,j))); if(nrm>0.0_dp)u(:,j)=u(:,j)/nrm; end do
    end if
    dr%d=d; dr%numdir=nd; dr%lambda=lam
    allocate(dr%evalues(size(ev)),dr%raw_evectors(d,size(ev)),dr%basis(d,nd),dr%direction(size(x,1),nd))
    allocate(dr%kernel(d,d),dr%kernel_mean(d,d),dr%kernel_cov(d,d))
    dr%evalues=ev; dr%raw_evectors=u; dr%basis=u(:,1:nd)
    dr%direction=matmul(x-spread(mu,1,size(x,1)),dr%basis)
    dr%kernel=m; dr%kernel_mean=mi; dr%kernel_cov=mii
    if(present(status)) status=0
  end subroutine fit_mclust_dr

  subroutine project_mclust_dr(dr,x,center,z)
    type(mclust_dr_fit),intent(in)::dr
    real(dp),intent(in)::x(:,:)
    real(dp),intent(in),optional::center(:)
    real(dp),allocatable,intent(out)::z(:,:)
    allocate(z(size(x,1),dr%numdir))
    if(present(center)) then; z=matmul(x-spread(center,1,size(x,1)),dr%basis); else; z=matmul(x,dr%basis); end if
  end subroutine project_mclust_dr

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:); real(dp)::c(size(a),size(b)); integer::i
    do i=1,size(a); c(i,:)=a(i)*b; end do
  end function outer
end module mclust_dr
