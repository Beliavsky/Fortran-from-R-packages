! SPDX-License-Identifier: MIT
module bekks_linalg
  use bekks_kinds, only: dp
  implicit none
  private
  public :: covariance_matrix, outer_product, kron, vec_col, mat_col
  public :: symmetric_inverse, general_inverse, solve_linear, logdet_quadratic
  public :: spectral_radius, symmetric_sqrt, cholesky_lower, matrix_power

  interface
    subroutine dpotrf(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine
    subroutine dpotri(uplo,n,a,lda,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: info
    end subroutine
    subroutine dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
      import dp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n,nrhs,lda,ldb
      real(dp), intent(in) :: a(lda,*)
      real(dp), intent(inout) :: b(ldb,*)
      integer, intent(out) :: info
    end subroutine
    subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      integer, intent(in) :: n,nrhs,lda,ldb
      integer, intent(out) :: ipiv(*)
      real(dp), intent(inout) :: a(lda,*),b(ldb,*)
      integer, intent(out) :: info
    end subroutine
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: w(*)
      integer, intent(out) :: info
    end subroutine
    subroutine dgeev(jobvl,jobvr,n,a,lda,wr,wi,vl,ldvl,vr,ldvr,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobvl,jobvr
      integer, intent(in) :: n,lda,ldvl,ldvr,lwork
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: wr(*),wi(*),vl(ldvl,*),vr(ldvr,*)
      integer, intent(out) :: info
    end subroutine
    subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobu,jobvt
      integer, intent(in) :: m,n,lda,ldu,ldvt,lwork
      real(dp), intent(inout) :: a(lda,*),work(*)
      real(dp), intent(out) :: s(*),u(ldu,*),vt(ldvt,*)
      integer, intent(out) :: info
    end subroutine
  end interface

contains

  function covariance_matrix(x) result(c)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable :: c(:,:)
    integer :: t
    t=size(x,1); allocate(c(size(x,2),size(x,2)))
    c=matmul(transpose(x),x)/real(max(1,t),dp)
  end function covariance_matrix

  pure function outer_product(x) result(a)
    real(dp), intent(in) :: x(:)
    real(dp) :: a(size(x),size(x))
    a=spread(x,2,size(x))*spread(x,1,size(x))
  end function outer_product

  pure function kron(a,b) result(c)
    real(dp), intent(in) :: a(:,:),b(:,:)
    real(dp) :: c(size(a,1)*size(b,1),size(a,2)*size(b,2))
    integer :: i,j,m,n
    m=size(b,1); n=size(b,2)
    do j=1,size(a,2)
      do i=1,size(a,1)
        c((i-1)*m+1:i*m,(j-1)*n+1:j*n)=a(i,j)*b
      end do
    end do
  end function kron

  pure function vec_col(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: v(size(a))
    v=reshape(a,[size(a)])
  end function vec_col

  pure function mat_col(v,n,m) result(a)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: n,m
    real(dp) :: a(n,m)
    a=reshape(v,[n,m])
  end function mat_col

  subroutine cholesky_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n,i,j
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

  subroutine symmetric_inverse(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n,i,j
    n=size(a,1); ainv=a
    call dpotrf('L',n,ainv,n,info)
    if(info/=0) return
    call dpotri('L',n,ainv,n,info)
    if(info/=0) return
    do j=1,n
      do i=1,j-1
        ainv(i,j)=ainv(j,i)
      end do
    end do
  end subroutine symmetric_inverse

  subroutine general_inverse(a,ainv,info,tol)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tol
    integer :: m,n,k,lwork,i
    real(dp), allocatable :: ac(:,:),s(:),u(:,:),vt(:,:),work(:)
    real(dp) :: threshold
    m=size(a,1); n=size(a,2); k=min(m,n)
    allocate(ac(m,n),s(k),u(m,m),vt(n,n),work(1)); ac=a
    lwork=-1
    call dgesvd('A','A',m,n,ac,m,s,u,m,vt,n,work,lwork,info)
    if(info/=0) return
    lwork=max(1,int(work(1))); deallocate(work); allocate(work(lwork)); ac=a
    call dgesvd('A','A',m,n,ac,m,s,u,m,vt,n,work,lwork,info)
    if(info/=0) return
    threshold=max(m,n)*epsilon(1.0_dp)*maxval(s)
    if(present(tol))threshold=tol
    ainv=0.0_dp
    do i=1,k
      if(s(i)>threshold) ainv=ainv+spread(vt(i,:),2,m)*spread(u(:,i),1,n)/s(i)
    end do
  end subroutine general_inverse

  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ac(:,:),bc(:,:)
    integer, allocatable :: ipiv(:)
    integer :: n
    n=size(a,1); allocate(ac(n,n),bc(n,1),ipiv(n)); ac=a; bc(:,1)=b
    call dgesv(n,1,ac,n,ipiv,bc,n,info)
    x=bc(:,1)
  end subroutine solve_linear

  subroutine logdet_quadratic(a,x,logdet,quad,info)
    real(dp), intent(in) :: a(:,:),x(:)
    real(dp), intent(out) :: logdet,quad
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:),b(:,:)
    integer :: n,i
    n=size(a,1); allocate(l(n,n),b(n,1)); l=0.5_dp*(a+transpose(a)); b(:,1)=x
    call dpotrf('L',n,l,n,info)
    if(info/=0) then
      logdet=huge(1.0_dp); quad=huge(1.0_dp); return
    end if
    logdet=0.0_dp
    do i=1,n
      if(l(i,i)<=0.0_dp) then
        info=1; return
      end if
      logdet=logdet+2.0_dp*log(l(i,i))
    end do
    call dpotrs('L',n,1,l,n,b,n,info)
    if(info==0) quad=dot_product(x,b(:,1))
  end subroutine logdet_quadratic

  real(dp) function spectral_radius(a,info) result(rho)
    real(dp), intent(in) :: a(:,:)
    integer, intent(out) :: info
    integer :: n,lwork
    real(dp), allocatable :: ac(:,:),wr(:),wi(:),work(:),vl(:,:),vr(:,:)
    n=size(a,1); allocate(ac(n,n),wr(n),wi(n),vl(1,1),vr(1,1),work(1)); ac=a
    lwork=-1
    call dgeev('N','N',n,ac,n,wr,wi,vl,1,vr,1,work,lwork,info)
    if(info/=0) then; rho=huge(1.0_dp); return; end if
    lwork=max(1,int(work(1))); deallocate(work); allocate(work(lwork)); ac=a
    call dgeev('N','N',n,ac,n,wr,wi,vl,1,vr,1,work,lwork,info)
    if(info==0) then
      rho=maxval(sqrt(wr*wr+wi*wi))
    else
      rho=huge(1.0_dp)
    end if
  end function spectral_radius

  subroutine symmetric_sqrt(a,sqrt_a,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: sqrt_a(:,:)
    integer, intent(out) :: info
    integer :: n,lwork,i
    real(dp), allocatable :: ac(:,:),w(:),work(:)
    n=size(a,1); allocate(ac(n,n),w(n),work(1)); ac=0.5_dp*(a+transpose(a))
    lwork=-1; call dsyev('V','L',n,ac,n,w,work,lwork,info)
    if(info/=0)return
    lwork=max(1,int(work(1))); deallocate(work); allocate(work(lwork)); ac=0.5_dp*(a+transpose(a))
    call dsyev('V','L',n,ac,n,w,work,lwork,info)
    if(info/=0)return
    sqrt_a=0.0_dp
    do i=1,n
      if(w(i)<-1.0e-10_dp) then
        info=1; return
      end if
      sqrt_a=sqrt_a+sqrt(max(w(i),0.0_dp))*outer_product(ac(:,i))
    end do
  end subroutine symmetric_sqrt

  function matrix_power(a,k) result(p)
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: k
    real(dp) :: p(size(a,1),size(a,2)),base(size(a,1),size(a,2))
    integer :: n,e,i
    n=size(a,1); p=0.0_dp
    do i=1,n; p(i,i)=1.0_dp; end do
    if(k<=0)return
    base=a; e=k
    do while(e>0)
      if(mod(e,2)==1)p=matmul(p,base)
      e=e/2
      if(e>0)base=matmul(base,base)
    end do
  end function matrix_power

end module bekks_linalg
