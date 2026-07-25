! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module sv_linalg
  use sv_kinds, only : dp, log2pi
  use sv_rng, only : randn
  implicit none
  private
  public :: chol_lower, solve_spd, inverse_spd, logdet_spd, mvn_logpdf, mvn_draw
  public :: symmetric_eigen, covariance_matrix, correlation_from_covariance
  interface
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine dpotrf
    subroutine dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
      import dp
      character(1), intent(in) :: uplo
      integer, intent(in) :: n,nrhs,lda,ldb
      real(dp), intent(in) :: a(lda,*)
      real(dp), intent(inout) :: b(ldb,*)
      integer, intent(out) :: info
    end subroutine dpotrs
    subroutine dpotri(uplo,n,a,lda,info)
      import dp
      character(1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine dpotri
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(1), intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: w(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface
contains
  subroutine chol_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n,i,j
    n=size(a,1); allocate(l(n,n)); l=a
    call dpotrf('L',n,l,n,info)
    if (info==0) then
      do j=1,n
        do i=1,j-1
          l(i,j)=0.0_dp
        end do
      end do
    end if
  end subroutine chol_lower

  subroutine solve_spd(a,b,x,info)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:),bb(:,:)
    integer :: n
    n=size(a,1); call chol_lower(a,l,info)
    allocate(x(n)); x=0.0_dp
    if (info/=0) return
    allocate(bb(n,1)); bb(:,1)=b
    call dpotrs('L',n,1,l,n,bb,n,info)
    if (info==0) x=bb(:,1)
  end subroutine solve_spd

  subroutine inverse_spd(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n,i,j
    call chol_lower(a,ainv,info); if(info/=0)return
    n=size(a,1); call dpotri('L',n,ainv,n,info)
    if(info==0) then
      do j=1,n; do i=1,j-1; ainv(i,j)=ainv(j,i); end do; end do
    end if
  end subroutine inverse_spd

  real(dp) function logdet_spd(a,info) result(v)
    real(dp), intent(in) :: a(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:)
    integer :: i
    call chol_lower(a,l,info); v=-huge(1.0_dp)
    if(info==0) then
      v=0.0_dp; do i=1,size(a,1); v=v+2.0_dp*log(l(i,i)); end do
    end if
  end function logdet_spd

  real(dp) function mvn_logpdf(x,mean,cov) result(v)
    real(dp), intent(in) :: x(:),mean(:),cov(:,:)
    real(dp), allocatable :: sol(:)
    integer :: info
    real(dp) :: ld
    call solve_spd(cov,x-mean,sol,info)
    if(info/=0) then; v=-huge(1.0_dp); return; end if
    ld=logdet_spd(cov,info)
    if(info/=0) then; v=-huge(1.0_dp); return; end if
    v=-0.5_dp*(size(x)*log2pi+ld+dot_product(x-mean,sol))
  end function mvn_logpdf

  subroutine mvn_draw(mean,cov,x,info)
    real(dp), intent(in) :: mean(:),cov(:,:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:),z(:)
    integer :: i
    call chol_lower(cov,l,info); if(info/=0)return
    allocate(z(size(mean))); do i=1,size(z); z(i)=randn(); end do
    x=mean+matmul(l,z)
  end subroutine mvn_draw

  subroutine symmetric_eigen(a,values,vectors,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:),vectors(:,:)
    integer, intent(out) :: info
    integer :: n,lwork
    real(dp), allocatable :: work(:)
    n=size(a,1); allocate(vectors(n,n),values(n)); vectors=a
    lwork=max(1,3*n-1); allocate(work(lwork))
    call dsyev('V','U',n,vectors,n,values,work,lwork,info)
  end subroutine symmetric_eigen

  subroutine covariance_matrix(x,cov,demean)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: cov(:,:)
    logical, intent(in), optional :: demean
    logical :: dm
    real(dp), allocatable :: xc(:,:),means(:)
    integer :: n,j
    n=size(x,1); dm=.true.; if(present(demean)) dm=demean
    allocate(xc,source=x)
    if(dm) then
      allocate(means(size(x,2)))
      do j=1,size(x,2); means(j)=sum(x(:,j))/real(n,dp); xc(:,j)=xc(:,j)-means(j); end do
    end if
    cov=matmul(transpose(xc),xc)/real(max(1,n-1),dp)
  end subroutine covariance_matrix

  subroutine correlation_from_covariance(cov,cor)
    real(dp), intent(in) :: cov(:,:)
    real(dp), intent(out) :: cor(:,:)
    integer :: i,j,n
    n=size(cov,1)
    do j=1,n; do i=1,n
      cor(i,j)=cov(i,j)/sqrt(max(cov(i,i)*cov(j,j),tiny(1.0_dp)))
    end do; end do
    do i=1,n; cor(i,i)=1.0_dp; end do
  end subroutine correlation_from_covariance
end module sv_linalg
