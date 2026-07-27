! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_linalg
  use fportfolio_kinds, only: dp
  implicit none
  private
  public :: solve_linear, inverse_matrix, symmetric_eigen, make_positive_definite, &
            quadratic_form, matrix_rank, cholesky_factor

  interface
    subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      integer, intent(in) :: n,nrhs,lda,ldb
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*), b(ldb,*)
    end subroutine dgesv
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: w(*)
    end subroutine dsyev
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
    end subroutine dpotrf
  end interface
contains
  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), bb(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n
    n=size(b)
    allocate(aa(n,n),bb(n,1),ipiv(n),x(n))
    aa=a; bb(:,1)=b
    call dgesv(n,1,aa,n,ipiv,bb,n,info)
    if (info == 0) then
      x=bb(:,1)
    else
      x=0.0_dp
    end if
  end subroutine solve_linear

  subroutine inverse_matrix(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), b(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n,i
    n=size(a,1)
    allocate(aa(n,n),b(n,n),ipiv(n),ainv(n,n))
    aa=a; b=0.0_dp
    do i=1,n; b(i,i)=1.0_dp; end do
    call dgesv(n,n,aa,n,ipiv,b,n,info)
    if (info == 0) then
      ainv=b
    else
      ainv=0.0_dp
    end if
  end subroutine inverse_matrix

  subroutine symmetric_eigen(a,values,vectors,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: work(:), tmp(:,:)
    integer :: n,lwork
    n=size(a,1); lwork=max(1,3*n-1)
    allocate(tmp(n,n),values(n),work(lwork),vectors(n,n))
    tmp=0.5_dp*(a+transpose(a))
    call dsyev('V','U',n,tmp,n,values,work,lwork,info)
    vectors=tmp
  end subroutine symmetric_eigen

  subroutine make_positive_definite(a,apd,floor_ratio,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: apd(:,:)
    real(dp), intent(in), optional :: floor_ratio
    integer, intent(out) :: info
    real(dp), allocatable :: vals(:),vecs(:,:)
    real(dp) :: fr, floorv
    integer :: i,n
    fr=1.0e-8_dp; if (present(floor_ratio)) fr=floor_ratio
    call symmetric_eigen(a,vals,vecs,info)
    n=size(a,1); allocate(apd(n,n)); apd=0.0_dp
    if (info /= 0) return
    floorv=max(maxval(abs(vals))*fr,1.0e-12_dp)
    vals=max(vals,floorv)
    do i=1,n
      apd=apd+vals(i)*spread(vecs(:,i),2,n)*spread(vecs(:,i),1,n)
    end do
    apd=0.5_dp*(apd+transpose(apd))
  end subroutine make_positive_definite

  pure real(dp) function quadratic_form(x,a) result(q)
    real(dp), intent(in) :: x(:),a(:,:)
    q=dot_product(x,matmul(a,x))
  end function quadratic_form

  subroutine matrix_rank(a,rank,tol)
    real(dp), intent(in) :: a(:,:)
    integer, intent(out) :: rank
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: ata(:,:),vals(:),vecs(:,:)
    real(dp) :: threshold
    integer :: info
    ata=matmul(transpose(a),a)
    call symmetric_eigen(ata,vals,vecs,info)
    if (info /= 0) then; rank=0; return; end if
    threshold=max(size(a,1),size(a,2))*epsilon(1.0_dp)*sqrt(maxval(max(vals,0.0_dp)))
    if (present(tol)) threshold=tol
    rank=count(sqrt(max(vals,0.0_dp))>threshold)
  end subroutine matrix_rank

  subroutine cholesky_factor(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n,i,j
    n=size(a,1); allocate(l(n,n)); l=0.5_dp*(a+transpose(a))
    call dpotrf('L',n,l,n,info)
    if (info == 0) then
      do j=1,n; do i=1,j-1; l(i,j)=0.0_dp; end do; end do
    end if
  end subroutine cholesky_factor
end module fportfolio_linalg
