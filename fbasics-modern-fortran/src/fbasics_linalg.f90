! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_linalg
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use fbasics_kinds, only: dp
  implicit none
  private
  public :: matrix_inverse, matrix_rank, matrix_norm, matrix_trace
  public :: hilbert_matrix, pascal_matrix, kronecker_product
  public :: lower_triangle, upper_triangle, vec_matrix, vech_matrix
  public :: is_positive_definite, make_positive_definite, grid_vector
  public :: lag_matrix, polynomial_distributed_lags
  interface
    subroutine dgetrf(m,n,a,lda,ipiv,info)
      import dp
      integer, intent(in) :: m,n,lda
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine
    subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
      import dp
      integer, intent(in) :: n,lda,lwork
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(in) :: ipiv(*)
      real(dp), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine
    subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobu,jobvt
      integer, intent(in) :: m,n,lda,ldu,ldvt,lwork
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: s(*),u(ldu,*),vt(ldvt,*),work(*)
      integer, intent(out) :: info
    end subroutine
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character(len=1), intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: w(*),work(*)
      integer, intent(out) :: info
    end subroutine
  end interface
contains
  subroutine matrix_inverse(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n,lwork
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: work(:)
    n=size(a,1)
    if (size(a,2)/=n) then; info=-1; allocate(ainv(0,0)); return; end if
    allocate(ainv(n,n),ipiv(n)); ainv=a
    call dgetrf(n,n,ainv,n,ipiv,info); if (info/=0) return
    lwork=max(1,64*n); allocate(work(lwork))
    call dgetri(n,ainv,n,ipiv,work,lwork,info)
  end subroutine matrix_inverse

  integer function matrix_rank(a, tol) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    integer :: m,n,mn,lwork,info
    real(dp) :: threshold
    real(dp), allocatable :: ac(:,:),s(:),u(:,:),vt(:,:),work(:)
    m=size(a,1); n=size(a,2); mn=min(m,n); lwork=max(1,8*max(m,n))
    allocate(ac(m,n),s(mn),u(1,1),vt(1,1),work(lwork)); ac=a
    call dgesvd('N','N',m,n,ac,m,s,u,1,vt,1,work,lwork,info)
    if (info/=0) then; r=0; return; end if
    threshold=max(m,n)*epsilon(1.0_dp)*maxval(s)
    if (present(tol)) threshold=tol
    r=count(s>threshold)
  end function matrix_rank

  real(dp) function matrix_norm(a,p) result(v)
    real(dp), intent(in) :: a(:,:)
    integer, intent(in) :: p
    integer :: m,n,mn,lwork,info
    real(dp), allocatable :: ac(:,:),s(:),u(:,:),vt(:,:),work(:)
    select case(p)
    case(1); v=maxval(sum(abs(a),dim=1))
    case(2)
      m=size(a,1); n=size(a,2); mn=min(m,n); lwork=max(1,8*max(m,n))
      allocate(ac(m,n),s(mn),u(1,1),vt(1,1),work(lwork)); ac=a
      call dgesvd('N','N',m,n,ac,m,s,u,1,vt,1,work,lwork,info)
      if (info==0) then; v=maxval(s); else; v=huge(1.0_dp); end if
    case default; v=maxval(sum(abs(a),dim=2))
    end select
  end function matrix_norm

  pure real(dp) function matrix_trace(a) result(v)
    real(dp), intent(in) :: a(:,:)
    integer :: i
    v=0.0_dp
    do i=1,min(size(a,1),size(a,2)); v=v+a(i,i); end do
  end function matrix_trace

  function hilbert_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp), allocatable :: a(:,:)
    integer :: i,j
    allocate(a(n,n)); do j=1,n; do i=1,n; a(i,j)=1.0_dp/real(i+j-1,dp); end do; end do
  end function hilbert_matrix

  function pascal_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp), allocatable :: a(:,:)
    integer :: i,j
    allocate(a(n,n)); do j=1,n; do i=1,n; a(i,j)=binom(i+j-2,j-1); end do; end do
  contains
    pure real(dp) function binom(nn,kk) result(v)
      integer,intent(in)::nn,kk
      integer::q,k0
      v=1.0_dp; k0=min(kk,nn-kk)
      do q=1,k0; v=v*real(nn-k0+q,dp)/real(q,dp); end do
    end function
  end function pascal_matrix

  function kronecker_product(a,b) result(k)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp),allocatable::k(:,:)
    integer::i,j,mb,nb
    mb=size(b,1); nb=size(b,2); allocate(k(size(a,1)*mb,size(a,2)*nb))
    do j=1,size(a,2); do i=1,size(a,1)
      k((i-1)*mb+1:i*mb,(j-1)*nb+1:j*nb)=a(i,j)*b
    end do; end do
  end function

  function lower_triangle(a) result(b)
    real(dp),intent(in)::a(:,:); real(dp),allocatable::b(:,:); integer::i,j
    allocate(b(size(a,1),size(a,2))); b=a
    do j=1,size(b,2); do i=1,min(j-1,size(b,1)); b(i,j)=0.0_dp; end do; end do
  end function
  function upper_triangle(a) result(b)
    real(dp),intent(in)::a(:,:); real(dp),allocatable::b(:,:); integer::i,j
    allocate(b(size(a,1),size(a,2))); b=a
    do j=1,size(b,2); do i=j+1,size(b,1); b(i,j)=0.0_dp; end do; end do
  end function
  function vec_matrix(a) result(v)
    real(dp),intent(in)::a(:,:); real(dp),allocatable::v(:)
    allocate(v(size(a))); v=reshape(a,[size(a)])
  end function
  function vech_matrix(a) result(v)
    real(dp),intent(in)::a(:,:); real(dp),allocatable::v(:); integer::i,j,k,n
    n=min(size(a,1),size(a,2)); allocate(v(n*(n+1)/2)); k=0
    do i=1,n; do j=1,i; k=k+1; v(k)=a(i,j); end do; end do
  end function

  logical function is_positive_definite(a,tol) result(ok)
    real(dp),intent(in)::a(:,:); real(dp),intent(in),optional::tol
    real(dp),allocatable::ac(:,:),w(:),work(:); real(dp) :: threshold
    integer::n,lwork,info
    n=size(a,1); if(size(a,2)/=n) then; ok=.false.; return; end if
    lwork=max(1,4*n); allocate(ac(n,n),w(n),work(lwork)); ac=0.5_dp*(a+transpose(a))
    call dsyev('N','U',n,ac,n,w,work,lwork,info)
    threshold=n*max(1.0_dp,maxval(abs(w)))*epsilon(1.0_dp); if(present(tol))threshold=tol
    ok=info==0 .and. minval(w)>threshold
  end function

  subroutine make_positive_definite(a,out,tol,info)
    real(dp),intent(in)::a(:,:); real(dp),allocatable,intent(out)::out(:,:)
    real(dp),intent(in),optional::tol; integer,intent(out)::info
    real(dp),allocatable::eigvec(:,:),w(:),work(:); real(dp) :: threshold
    integer::n,lwork,i
    n=size(a,1); if(size(a,2)/=n) then; info=-1; allocate(out(0,0)); return; end if
    lwork=max(1,4*n); allocate(eigvec(n,n),w(n),work(lwork)); eigvec=0.5_dp*(a+transpose(a))
    call dsyev('V','U',n,eigvec,n,w,work,lwork,info); if(info/=0) then; allocate(out(0,0)); return; end if
    threshold=2.0_dp*n*max(1.0_dp,maxval(abs(w)))*epsilon(1.0_dp); if(present(tol))threshold=2.0_dp*tol
    w=max(w,threshold); allocate(out(n,n)); out=0.0_dp
    do i=1,n; out=out+w(i)*outer(eigvec(:,i),eigvec(:,i)); end do
  contains
    pure function outer(x,y) result(z)
      real(dp),intent(in)::x(:),y(:); real(dp)::z(size(x),size(y)); integer::j
      do j=1,size(y); z(:,j)=x*y(j); end do
    end function
  end subroutine

  subroutine grid_vector(x,y,gx,gy)
    real(dp),intent(in)::x(:),y(:); real(dp),allocatable,intent(out)::gx(:),gy(:)
    integer::i,j,k
    allocate(gx(size(x)*size(y)),gy(size(x)*size(y))); k=0
    do j=1,size(y); do i=1,size(x); k=k+1; gx(k)=x(i); gy(k)=y(j); end do; end do
  end subroutine

  function lag_matrix(x,lags,trim) result(m)
    real(dp),intent(in)::x(:); integer,intent(in)::lags(:); logical,intent(in),optional::trim
    real(dp),allocatable::m(:,:); real(dp) :: nanv; integer::i,j,k,n,first,last
    logical::do_trim
    n=size(x); nanv=ieee_value(0.0_dp,ieee_quiet_nan); do_trim=.false.; if(present(trim))do_trim=trim
    allocate(m(n,size(lags))); m=nanv
    do j=1,size(lags)
      k=lags(j)
      do i=1,n
        if(i-k>=1 .and. i-k<=n)m(i,j)=x(i-k)
      end do
    end do
    if(do_trim) then
      first=1+max(0,maxval(lags)); last=n+min(0,minval(lags))
      if(last>=first)m=m(first:last,:)
    end if
  end function

  function polynomial_distributed_lags(x,d,q,trim) result(z)
    real(dp),intent(in)::x(:); integer,intent(in)::d,q; logical,intent(in),optional::trim
    real(dp),allocatable::z(:,:),m(:,:); integer,allocatable::lags(:); integer::i,j,k,n,first
    logical::do_trim
    if(q<=d)then; allocate(z(0,0)); return; end if
    allocate(lags(q)); lags=[(i,i=1,q)]; m=lag_matrix(x,lags,.false.); n=size(x)
    allocate(z(n,d+1)); z=0.0_dp
    do i=1,n
      do j=0,d
        if(i<=q) then
          z(i,j+1)=ieee_value(0.0_dp,ieee_quiet_nan)
        else
          z(i,j+1)=sum([(real(k,dp)**j*m(i,k),k=1,q)])
          if(j==0)z(i,j+1)=z(i,j+1)+x(i)
        end if
      end do
    end do
    do_trim=.false.; if(present(trim))do_trim=trim
    if(do_trim) then; first=q+1; z=z(first:n,:); end if
  end function
end module fbasics_linalg
