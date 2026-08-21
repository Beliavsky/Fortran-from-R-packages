! SPDX-License-Identifier: GPL-2.0-or-later
module compositions_linalg
  use compositions_kinds, only: dp
  implicit none
  private
  public :: symmetric_eigen, pseudoinverse, invert_matrix, determinant_spd
  public :: covariance_matrix, solve_least_squares, matrix_rank, chol_lower

  interface
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: w(*)
      integer, intent(out) :: info
    end subroutine dsyev
    subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobu,jobvt
      integer, intent(in) :: m,n,lda,ldu,ldvt,lwork
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: s(*),u(ldu,*),vt(ldvt,*)
      integer, intent(out) :: info
    end subroutine dgesvd
    subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      integer, intent(in) :: n,nrhs,lda,ldb
      real(dp), intent(inout) :: a(lda,*),b(ldb,*)
      integer, intent(out) :: ipiv(*),info
    end subroutine dgesv
    subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
      import dp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m,n,nrhs,lda,ldb,lwork
      real(dp), intent(inout) :: a(lda,*),b(ldb,*),work(*)
      integer, intent(out) :: info
    end subroutine dgels
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine dpotrf
  end interface
contains
  subroutine symmetric_eigen(a, values, vectors, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: tmp(:,:), work(:)
    integer :: n, lwork, ierr
    n = size(a,1)
    if(size(a,2) /= n) error stop 'symmetric_eigen: matrix must be square'
    allocate(tmp(n,n), values(n))
    tmp = a
    lwork = max(1, 3*n-1)
    allocate(work(lwork))
    call dsyev('V','U',n,tmp,n,values,work,lwork,ierr)
    if(present(info)) info = ierr
    if(ierr /= 0) error stop 'symmetric_eigen: dsyev failed'
    vectors = tmp
  end subroutine symmetric_eigen

  subroutine pseudoinverse(a, ainv, tol)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: ac(:,:), s(:), u(:,:), vt(:,:), work(:), sinv(:)
    real(dp) :: threshold
    integer :: m,n,k,lwork,info,i
    m=size(a,1); n=size(a,2); k=min(m,n)
    allocate(ac(m,n),s(k),u(m,m),vt(n,n))
    ac=a
    lwork=max(1,5*max(m,n))
    allocate(work(lwork))
    call dgesvd('A','A',m,n,ac,m,s,u,m,vt,n,work,lwork,info)
    if(info/=0) error stop 'pseudoinverse: dgesvd failed'
    threshold = 1.0e-10_dp
    if(present(tol)) threshold=tol
    if(k>0) threshold=threshold*maxval(s)
    allocate(sinv(k)); sinv=0.0_dp
    do i=1,k
      if(s(i)>threshold) sinv(i)=1.0_dp/s(i)
    end do
    allocate(ainv(n,m)); ainv=0.0_dp
    do i=1,k
      if(sinv(i)/=0.0_dp) ainv = ainv + sinv(i)*outer_product(vt(i,:),u(:,i))
    end do
  contains
    function outer_product(x,y) result(z)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: z(size(x),size(y))
      integer :: ii,jj
      do jj=1,size(y); do ii=1,size(x); z(ii,jj)=x(ii)*y(jj); end do; end do
    end function outer_product
  end subroutine pseudoinverse

  subroutine invert_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: ac(:,:), b(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n,i,ierr
    n=size(a,1)
    if(size(a,2)/=n) error stop 'invert_matrix: square matrix required'
    allocate(ac(n,n),b(n,n),ipiv(n)); ac=a; b=0.0_dp
    do i=1,n; b(i,i)=1.0_dp; end do
    call dgesv(n,n,ac,n,ipiv,b,n,ierr)
    if(present(info)) info=ierr
    if(ierr/=0) error stop 'invert_matrix: singular matrix'
    ainv=b
  end subroutine invert_matrix

  function determinant_spd(a) result(det)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: det
    real(dp), allocatable :: c(:,:)
    integer :: n,i,info
    n=size(a,1)
    if(size(a,2)/=n) error stop 'determinant_spd: square matrix required'
    allocate(c(n,n)); c=a
    call dpotrf('L',n,c,n,info)
    if(info/=0) then
      det=0.0_dp
      return
    end if
    det=1.0_dp
    do i=1,n; det=det*c(i,i)*c(i,i); end do
  end function determinant_spd

  subroutine chol_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out), optional :: info
    integer :: n,i,j,ierr
    n=size(a,1)
    allocate(l(n,n)); l=a
    call dpotrf('L',n,l,n,ierr)
    do j=1,n; do i=1,n; if(i<j) l(i,j)=0.0_dp; end do; end do
    if(present(info)) info=ierr
    if(ierr/=0) error stop 'chol_lower: matrix is not positive definite'
  end subroutine chol_lower

  subroutine covariance_matrix(x,cov,center)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    real(dp), allocatable, intent(out), optional :: center(:)
    real(dp), allocatable :: mu(:),xc(:,:)
    integer :: n,p
    n=size(x,1); p=size(x,2)
    allocate(mu(p)); mu=sum(x,dim=1)/real(n,dp)
    allocate(xc(n,p)); xc=x-spread(mu,1,n)
    allocate(cov(p,p))
    if(n>1) then
      cov=matmul(transpose(xc),xc)/real(n-1,dp)
    else
      cov=0.0_dp
    end if
    if(present(center)) center=mu
  end subroutine covariance_matrix

  subroutine solve_least_squares(x,y,beta,residuals,info)
    real(dp), intent(in) :: x(:,:),y(:,:)
    real(dp), allocatable, intent(out) :: beta(:,:),residuals(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: a(:,:),b(:,:),work(:)
    integer :: m,n,nrhs,ldb,lwork,ierr
    m=size(x,1); n=size(x,2); nrhs=size(y,2); ldb=max(m,n)
    if(size(y,1)/=m) error stop 'solve_least_squares: row mismatch'
    allocate(a(m,n),b(ldb,nrhs)); a=x; b=0.0_dp; b(1:m,:)=y
    lwork=max(1,4*max(m,n)+max(m,nrhs))
    allocate(work(lwork))
    call dgels('N',m,n,nrhs,a,m,b,ldb,work,lwork,ierr)
    if(present(info)) info=ierr
    if(ierr/=0) error stop 'solve_least_squares: dgels failed'
    beta=b(1:n,:)
    residuals=y-matmul(x,beta)
  end subroutine solve_least_squares

  integer function matrix_rank(a,tol) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: pinv(:,:),ac(:,:),s(:),u(:,:),vt(:,:),work(:)
    real(dp) :: th
    integer :: m,n,k,lwork,info
    m=size(a,1); n=size(a,2); k=min(m,n)
    allocate(ac(m,n),s(k),u(1,1),vt(1,1)); ac=a
    lwork=max(1,5*max(m,n)); allocate(work(lwork))
    call dgesvd('N','N',m,n,ac,m,s,u,1,vt,1,work,lwork,info)
    if(info/=0) error stop 'matrix_rank: dgesvd failed'
    th=1.0e-10_dp
    if(present(tol)) th=tol
    if(k>0) th=th*maxval(s)
    r=count(s>th)
  end function matrix_rank
end module compositions_linalg
