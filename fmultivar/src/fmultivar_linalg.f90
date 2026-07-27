! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_linalg
  use fmultivar_kinds, only : dp
  implicit none
  private
  public :: cholesky_lower, inverse_spd, inverse_general, logdet_spd, solve_spd, sample_mean_cov
  public :: is_symmetric, is_positive_definite, covariance_to_correlation

  interface
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine dpotrf
    subroutine dpotri(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine dpotri
    subroutine dgetrf(m,n,a,lda,ipiv,info)
      import dp
      integer, intent(in) :: m,n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine dgetrf
    subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      import dp
      integer, intent(in) :: n,lda,lwork
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(in) :: ipiv(*)
      real(dp), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine dgetri
    subroutine dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,nrhs,lda,ldb
      real(dp), intent(in) :: a(lda,*)
      real(dp), intent(inout) :: b(ldb,*)
      integer, intent(out) :: info
    end subroutine dpotrs
  end interface
contains
  subroutine cholesky_lower(a, l, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: n, info, i
    n = size(a,1)
    allocate(l(n,n)); l = a
    call dpotrf('L',n,l,n,info)
    ok = info == 0
    if (.not.ok) then; l = 0.0_dp; return; end if
    do i=1,n-1; l(i,i+1:n)=0.0_dp; end do
  end subroutine cholesky_lower

  subroutine inverse_spd(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    integer :: n, info, i, j
    n=size(a,1); allocate(ainv(n,n)); ainv=a
    call dpotrf('L',n,ainv,n,info)
    if (info/=0) then; ok=.false.; ainv=0.0_dp; return; end if
    call dpotri('L',n,ainv,n,info)
    ok=info==0
    if (.not.ok) then; ainv=0.0_dp; return; end if
    do i=1,n; do j=i+1,n; ainv(i,j)=ainv(j,i); end do; end do
  end subroutine inverse_spd


  subroutine inverse_general(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    integer :: n, info, lwork
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: work(:)
    n = size(a,1)
    allocate(ainv(n,n),ipiv(n))
    ainv = a
    call dgetrf(n,n,ainv,n,ipiv,info)
    if (info /= 0) then
      ok = .false.; ainv = 0.0_dp; return
    end if
    lwork = max(1,64*n)
    allocate(work(lwork))
    call dgetri(n,ainv,n,ipiv,work,lwork,info)
    ok = info == 0
    if (.not.ok) ainv = 0.0_dp
  end subroutine inverse_general

  function logdet_spd(a, ok) result(v)
    real(dp), intent(in) :: a(:,:)
    logical, intent(out) :: ok
    real(dp) :: v
    real(dp), allocatable :: l(:,:)
    integer :: i
    call cholesky_lower(a,l,ok)
    if (.not.ok) then; v=huge(1.0_dp); return; end if
    v=0.0_dp
    do i=1,size(l,1); v=v+2.0_dp*log(l(i,i)); end do
  end function logdet_spd

  subroutine solve_spd(a,b,x,ok)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: l(:,:)
    integer :: n, nrhs, info
    n=size(a,1); nrhs=size(b,2)
    call cholesky_lower(a,l,ok)
    allocate(x(size(b,1),nrhs)); x=b
    if (.not.ok) then; x=0.0_dp; return; end if
    call dpotrs('L',n,nrhs,l,n,x,n,info)
    ok=info==0
  end subroutine solve_spd

  subroutine sample_mean_cov(x, mean, cov)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mean(:), cov(:,:)
    integer :: n,p,i
    n=size(x,1); p=size(x,2)
    allocate(mean(p),cov(p,p))
    mean=sum(x,dim=1)/real(n,dp)
    cov=0.0_dp
    do i=1,n
      cov=cov+matmul(reshape(x(i,:)-mean,[p,1]),reshape(x(i,:)-mean,[1,p]))
    end do
    cov=cov/real(max(1,n-1),dp)
  end subroutine sample_mean_cov

  pure function is_symmetric(a,tol) result(ok)
    real(dp), intent(in) :: a(:,:), tol
    logical :: ok
    ok=size(a,1)==size(a,2)
    if (ok) ok=maxval(abs(a-transpose(a)))<=tol
  end function is_symmetric

  function is_positive_definite(a) result(ok)
    real(dp), intent(in) :: a(:,:)
    logical :: ok
    real(dp), allocatable :: l(:,:)
    call cholesky_lower(a,l,ok)
  end function is_positive_definite

  subroutine covariance_to_correlation(cov, cor, scale)
    real(dp), intent(in) :: cov(:,:)
    real(dp), allocatable, intent(out) :: cor(:,:), scale(:)
    integer :: p,i,j
    p=size(cov,1); allocate(cor(p,p),scale(p))
    do i=1,p; scale(i)=sqrt(max(cov(i,i),tiny(1.0_dp))); end do
    do i=1,p; do j=1,p; cor(i,j)=cov(i,j)/(scale(i)*scale(j)); end do; end do
  end subroutine covariance_to_correlation
end module fmultivar_linalg
