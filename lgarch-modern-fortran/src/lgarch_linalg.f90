! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_linalg
  use lgarch_kinds, only : dp
  implicit none
  private
  public :: solve_linear, inverse_matrix, cholesky_lower, logdet_spd
  public :: symmetric_eigenvalues, spectral_radius, covariance_matrix, correlation_matrix
  interface
    subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      integer, intent(in) :: n,nrhs,lda,ldb
      integer, intent(out) :: ipiv(*),info
      real(dp), intent(inout) :: a(lda,*),b(ldb,*)
    end subroutine dgesv
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
    end subroutine dpotrf
    subroutine dgetrf(m,n,a,lda,ipiv,info)
      import dp
      integer, intent(in) :: m,n,lda
      integer, intent(out) :: ipiv(*),info
      real(dp), intent(inout) :: a(lda,*)
    end subroutine dgetrf
    subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      import dp
      integer, intent(in) :: n,lda,lwork
      integer, intent(in) :: ipiv(*)
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*),work(*)
    end subroutine dgetri
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: w(*)
    end subroutine dsyev
    subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobvl,jobvr
      integer, intent(in) :: n,lda,ldvl,ldvr,lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: wr(*),wi(*),vl(ldvl,*),vr(ldvr,*)
    end subroutine dgeev
  end interface
contains
  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ac(:,:), bc(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n
    n=size(a,1)
    if (size(a,2)/=n .or. size(b)/=n .or. size(x)/=n) error stop "solve_linear: size mismatch"
    allocate(ac(n,n),bc(n,1),ipiv(n)); ac=a; bc(:,1)=b
    call dgesv(n,1,ac,n,ipiv,bc,n,info)
    if (info==0) x=bc(:,1)
  end subroutine solve_linear

  subroutine inverse_matrix(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n,lwork
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: work(:)
    n=size(a,1)
    if(size(a,2)/=n .or. any(shape(ainv)/=shape(a))) error stop "inverse_matrix: size mismatch"
    ainv=a; allocate(ipiv(n)); call dgetrf(n,n,ainv,n,ipiv,info); if(info/=0)return
    lwork=max(1,64*n); allocate(work(lwork)); call dgetri(n,ainv,n,ipiv,work,lwork,info)
  end subroutine inverse_matrix

  subroutine cholesky_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: i,j,n
    n=size(a,1); l=a
    call dpotrf('L',n,l,n,info)
    if(info==0) then
      do j=1,n
        do i=1,j-1
          l(i,j)=0.0_dp
        end do
      end do
    end if
  end subroutine cholesky_lower

  subroutine logdet_spd(a,logdet,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:)
    integer :: i,n
    n=size(a,1); allocate(l(n,n)); call cholesky_lower(a,l,info)
    if(info/=0) then; logdet=huge(1.0_dp); return; end if
    logdet=0.0_dp
    do i=1,n; logdet=logdet+2.0_dp*log(l(i,i)); end do
  end subroutine logdet_spd

  subroutine symmetric_eigenvalues(a,w,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: w(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ac(:,:),work(:)
    integer :: n,lwork
    n=size(a,1); allocate(ac(n,n)); ac=a; lwork=max(1,3*n); allocate(work(lwork))
    call dsyev('N','U',n,ac,n,w,work,lwork,info)
  end subroutine symmetric_eigenvalues

  real(dp) function spectral_radius(phi) result(rho)
    real(dp), intent(in) :: phi(:)
    real(dp), allocatable :: a(:,:),wr(:),wi(:),vl(:,:),vr(:,:),work(:)
    integer :: p,lwork,info,i
    p=size(phi)
    if(p==0) then; rho=0.0_dp; return; end if
    allocate(a(p,p),wr(p),wi(p),vl(1,1),vr(1,1)); a=0.0_dp; a(1,:)=phi
    do i=2,p; a(i,i-1)=1.0_dp; end do
    lwork=max(1,4*p); allocate(work(lwork))
    call dgeev('N','N',p,a,p,wr,wi,vl,1,vr,1,work,lwork,info)
    if(info/=0) then; rho=huge(1.0_dp); else; rho=maxval(sqrt(wr*wr+wi*wi)); end if
  end function spectral_radius

  function covariance_matrix(x) result(cov)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: cov(:,:)
    real(dp), allocatable :: xc(:,:)
    real(dp) :: means(size(x,2))
    integer :: j,n
    n=size(x,1); allocate(cov(size(x,2),size(x,2)),xc(n,size(x,2)))
    do j=1,size(x,2); means(j)=sum(x(:,j))/real(n,dp); xc(:,j)=x(:,j)-means(j); end do
    if(n>1) then; cov=matmul(transpose(xc),xc)/real(n-1,dp); else; cov=0.0_dp; end if
  end function covariance_matrix

  function correlation_matrix(x) result(cor)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: cor(:,:),cov(:,:)
    real(dp) :: d(size(x,2))
    integer :: i,j
    cov=covariance_matrix(x); allocate(cor(size(cov,1),size(cov,2)))
    do i=1,size(d); d(i)=sqrt(max(cov(i,i),tiny(1.0_dp))); end do
    do j=1,size(d); do i=1,size(d); cor(i,j)=cov(i,j)/(d(i)*d(j)); end do; end do
  end function correlation_matrix
end module lgarch_linalg
