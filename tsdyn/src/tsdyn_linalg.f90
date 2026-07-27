! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_linalg
  use tsdyn_kinds, only: dp
  implicit none
  private
  public :: ols_fit, inverse_matrix, cholesky_lower, symmetric_eigen
  public :: general_eigenvalues, general_eigen, matrix_rank, covariance_matrix

  interface
    subroutine dgelsy(m,n,nrhs,a,lda,b,ldb,jpvt,rcond,rank,work,lwork,info)
      import dp
      integer, intent(in) :: m,n,nrhs,lda,ldb,lwork
      real(dp), intent(inout) :: a(lda,*), b(ldb,*)
      integer, intent(inout) :: jpvt(*)
      real(dp), intent(in) :: rcond
      integer, intent(out) :: rank,info
      real(dp), intent(inout) :: work(*)
    end subroutine dgelsy
    subroutine dgetrf(m,n,a,lda,ipiv,info)
      import dp
      integer,intent(in)::m,n,lda
      real(dp),intent(inout)::a(lda,*)
      integer,intent(out)::ipiv(*),info
    end subroutine dgetrf
    subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      import dp
      integer,intent(in)::n,lda,lwork
      real(dp),intent(inout)::a(lda,*),work(*)
      integer,intent(in)::ipiv(*)
      integer,intent(out)::info
    end subroutine dgetri
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(len=1),intent(in)::uplo
      integer,intent(in)::n,lda
      real(dp),intent(inout)::a(lda,*)
      integer,intent(out)::info
    end subroutine dpotrf
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(len=1),intent(in)::jobz,uplo
      integer,intent(in)::n,lda,lwork
      real(dp),intent(inout)::a(lda,*),work(*)
      real(dp),intent(out)::w(*)
      integer,intent(out)::info
    end subroutine dsyev
    subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
      import dp
      character(len=1),intent(in)::jobvl,jobvr
      integer,intent(in)::n,lda,ldvl,ldvr,lwork
      real(dp),intent(inout)::a(lda,*),work(*)
      real(dp),intent(out)::wr(*),wi(*),vl(ldvl,*),vr(ldvr,*)
      integer,intent(out)::info
    end subroutine dgeev
  end interface
contains
  subroutine ols_fit(x, y, beta, fitted, residuals, rank, ssr, info)
    real(dp), intent(in) :: x(:,:), y(:,:)
    real(dp), allocatable, intent(out) :: beta(:,:), fitted(:,:), residuals(:,:)
    integer, intent(out) :: rank, info
    real(dp), intent(out) :: ssr
    integer :: n, p, ny, ldb, lwork
    integer, allocatable :: jpvt(:)
    real(dp), allocatable :: a(:,:), b(:,:), work(:)
    real(dp) :: work_query(1)

    n = size(x,1); p = size(x,2); ny = size(y,2)
    if (size(y,1) /= n .or. n < 1 .or. p < 1) then
      info = -1; rank = 0; ssr = huge(1.0_dp)
      allocate(beta(0,0), fitted(0,0), residuals(0,0))
      return
    end if
    ldb = max(n,p)
    allocate(a(n,p), b(ldb,ny), jpvt(p))
    a = x; b = 0.0_dp; b(1:n,:) = y; jpvt = 0
    lwork = -1
    call dgelsy(n,p,ny,a,n,b,ldb,jpvt,1.0e-12_dp,rank,work_query,lwork,info)
    if (info /= 0) then
      allocate(beta(0,0), fitted(0,0), residuals(0,0))
      ssr = huge(1.0_dp)
      return
    end if
    lwork = max(1,int(work_query(1)))
    allocate(work(lwork))
    a = x; b = 0.0_dp; b(1:n,:) = y; jpvt = 0
    call dgelsy(n,p,ny,a,n,b,ldb,jpvt,1.0e-12_dp,rank,work,lwork,info)
    allocate(beta(p,ny), fitted(n,ny), residuals(n,ny))
    if (info /= 0) then
      beta = 0.0_dp; fitted = 0.0_dp; residuals = y
      ssr = huge(1.0_dp)
      return
    end if
    beta = b(1:p,:)
    fitted = matmul(x,beta)
    residuals = y - fitted
    ssr = sum(residuals*residuals)
  end subroutine ols_fit

  subroutine inverse_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n, lwork
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: work(:)
    real(dp) :: query(1)
    n = size(a,1)
    if (size(a,2) /= n .or. n < 1) then
      info = -1; allocate(ainv(0,0)); return
    end if
    allocate(ainv(n,n),ipiv(n))
    ainv = a
    call dgetrf(n,n,ainv,n,ipiv,info)
    if (info /= 0) return
    call dgetri(n,ainv,n,ipiv,query,-1,info)
    if (info /= 0) return
    lwork=max(1,int(query(1))); allocate(work(lwork))
    call dgetri(n,ainv,n,ipiv,work,lwork,info)
  end subroutine inverse_matrix

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n,i,j
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(l(0,0)); return
    end if
    allocate(l(n,n)); l=a
    call dpotrf('L',n,l,n,info)
    if(info==0)then
      do j=1,n
        do i=1,j-1
          l(i,j)=0.0_dp
        end do
      end do
    end if
  end subroutine cholesky_lower

  subroutine symmetric_eigen(a, values, vectors, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out) :: info
    integer :: n,lwork
    real(dp), allocatable :: work(:)
    real(dp) :: query(1)
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(values(0),vectors(0,0)); return
    end if
    allocate(values(n),vectors(n,n)); vectors=0.5_dp*(a+transpose(a))
    call dsyev('V','U',n,vectors,n,values,query,-1,info)
    if(info/=0)return
    lwork=max(1,int(query(1))); allocate(work(lwork))
    call dsyev('V','U',n,vectors,n,values,work,lwork,info)
  end subroutine symmetric_eigen

  subroutine general_eigenvalues(a, wr, wi, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: wr(:), wi(:)
    integer, intent(out) :: info
    integer :: n,lwork
    real(dp), allocatable :: ac(:,:),work(:),vl(:,:),vr(:,:)
    real(dp) :: query(1)
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(wr(0),wi(0)); return
    end if
    allocate(ac(n,n),wr(n),wi(n),vl(1,1),vr(1,1)); ac=a
    call dgeev('N','N',n,ac,n,wr,wi,vl,1,vr,1,query,-1,info)
    if(info/=0)return
    lwork=max(1,int(query(1))); allocate(work(lwork)); ac=a
    call dgeev('N','N',n,ac,n,wr,wi,vl,1,vr,1,work,lwork,info)
  end subroutine general_eigenvalues

  subroutine general_eigen(a, wr, wi, vr, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: wr(:), wi(:), vr(:,:)
    integer, intent(out) :: info
    integer :: n,lwork
    real(dp), allocatable :: ac(:,:),work(:),vl(:,:)
    real(dp) :: query(1)
    n=size(a,1)
    if(size(a,2)/=n)then
      info=-1; allocate(wr(0),wi(0),vr(0,0)); return
    end if
    allocate(ac(n,n),wr(n),wi(n),vl(1,1),vr(n,n)); ac=a
    call dgeev('N','V',n,ac,n,wr,wi,vl,1,vr,n,query,-1,info)
    if(info/=0)return
    lwork=max(1,int(query(1))); allocate(work(lwork)); ac=a
    call dgeev('N','V',n,ac,n,wr,wi,vl,1,vr,n,work,lwork,info)
  end subroutine general_eigen

  integer function matrix_rank(a, tol) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: beta(:,:),fit(:,:),res(:,:), y(:,:)
    real(dp)::ss
    integer::info
    allocate(y(size(a,1),1)); y=0.0_dp
    call ols_fit(a,y,beta,fit,res,r,ss,info)
    if(present(tol)) r = r + merge(0,0,tol>=0.0_dp)
  end function matrix_rank

  subroutine covariance_matrix(x, cov, center)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    logical, intent(in), optional :: center
    real(dp), allocatable :: xc(:,:),mu(:)
    logical :: do_center
    integer :: n,j
    n=size(x,1); do_center=.true.; if(present(center))do_center=center
    allocate(xc(size(x,1),size(x,2))); xc=x
    if(do_center)then
      allocate(mu(size(x,2))); mu=sum(x,dim=1)/real(max(1,n),dp)
      do j=1,size(x,2); xc(:,j)=xc(:,j)-mu(j); end do
    end if
    allocate(cov(size(x,2),size(x,2)))
    cov=matmul(transpose(xc),xc)/real(max(1,n-1),dp)
  end subroutine covariance_matrix
end module tsdyn_linalg
