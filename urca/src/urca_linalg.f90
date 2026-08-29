module urca_linalg
   use urca_kinds, only : dp
   implicit none
   private
   public :: invert_matrix, invert_spd, solve_matrix, solve_vector
   public :: symmetric_eigen, chol_lower, orthogonal_complement
   public :: eye, outer_product, trace_matrix, logdet_spd, determinant_spd

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
         integer, intent(out) :: ipiv(*), info
      end subroutine dgesv
      subroutine dgetrf(m, n, a, lda, ipiv, info)
         import dp
         integer, intent(in) :: m, n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: ipiv(*), info
      end subroutine dgetrf
      subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
         import dp
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*), work(*)
         integer, intent(in) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine dgetri
      subroutine dpotrf(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotrf
      subroutine dpotri(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotri
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*), work(*)
         real(dp), intent(out) :: w(*)
         integer, intent(out) :: info
      end subroutine dsyev
      subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
         import dp
         integer, intent(in) :: m, n, lda, lwork
         real(dp), intent(inout) :: a(lda,*), work(*)
         real(dp), intent(out) :: tau(*)
         integer, intent(out) :: info
      end subroutine dgeqrf
      subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
         import dp
         integer, intent(in) :: m, n, k, lda, lwork
         real(dp), intent(inout) :: a(lda,*), work(*)
         real(dp), intent(in) :: tau(*)
         integer, intent(out) :: info
      end subroutine dorgqr
   end interface
contains
   function eye(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function eye

   pure function outer_product(a,b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a),size(b))
      integer :: i
      do i = 1, size(a)
         c(i,:) = a(i)*b
      end do
   end function outer_product

   pure real(dp) function trace_matrix(a) result(v)
      real(dp), intent(in) :: a(:,:)
      integer :: i
      v = 0.0_dp
      do i=1,min(size(a,1),size(a,2))
         v=v+a(i,i)
      end do
   end function trace_matrix

   subroutine invert_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      integer :: n, lwork
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: work(:)
      real(dp) :: q(1)
      n=size(a,1)
      allocate(ainv(n,n),ipiv(n))
      if(size(a,2)/=n) then
      info=-1
      ainv=0
      return
      end if
      ainv=a
      call dgetrf(n,n,ainv,n,ipiv,info)
      if(info/=0)return
      call dgetri(n,ainv,n,ipiv,q,-1,info)
      if(info/=0)return
      lwork=max(1,int(q(1)))
      allocate(work(lwork))
      call dgetri(n,ainv,n,ipiv,work,lwork,info)
   end subroutine invert_matrix

   subroutine invert_spd(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      integer :: n,i,j
      n=size(a,1)
      allocate(ainv(n,n))
      ainv=a
      if(size(a,2)/=n) then
      info=-1
      return
      end if
      call dpotrf('L',n,ainv,n,info)
      if(info/=0)return
      call dpotri('L',n,ainv,n,info)
      if(info/=0)return
      do j=1,n
      do i=1,j-1
      ainv(i,j)=ainv(j,i)
      end do
      end do
   end subroutine invert_spd

   subroutine solve_matrix(a,b,x,info)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: ac(:,:),bc(:,:)
      integer, allocatable :: ipiv(:)
      integer :: n,nrhs
      n=size(a,1)
      nrhs=size(b,2)
      allocate(ac(n,n),bc(n,nrhs),ipiv(n),x(n,nrhs))
      if(size(a,2)/=n .or. size(b,1)/=n) then
      info=-1
      x=0
      return
      end if
      ac=a
      bc=b
      call dgesv(n,nrhs,ac,n,ipiv,bc,n,info)
      x=bc
   end subroutine solve_matrix

   subroutine solve_vector(a,b,x,info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: bm(:,:),xm(:,:)
      allocate(bm(size(b),1))
      bm(:,1)=b
      call solve_matrix(a,bm,xm,info)
      allocate(x(size(b)))
      if(size(xm,1)==size(b)) x=xm(:,1)
   end subroutine solve_vector

   subroutine chol_lower(a,l,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      integer, intent(out) :: info
      integer :: n,i,j
      n=size(a,1)
      allocate(l(n,n))
      l=a
      call dpotrf('L',n,l,n,info)
      if(info/=0)return
      do j=1,n
      do i=1,j-1
      l(i,j)=0.0_dp
      end do
      end do
   end subroutine chol_lower

   subroutine symmetric_eigen(a, values, vectors, info, descending)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:),vectors(:,:)
      integer, intent(out) :: info
      logical, intent(in), optional :: descending
      real(dp), allocatable :: ac(:,:),work(:),tmp(:)
      real(dp) :: q(1),tv
      integer :: n,lwork,i,j
      logical :: desc
      n=size(a,1)
      allocate(ac(n,n),values(n),vectors(n,n))
      ac=0.5_dp*(a+transpose(a))
      call dsyev('V','U',n,ac,n,values,q,-1,info)
      if(info/=0)return
      lwork=max(1,int(q(1)))
      allocate(work(lwork))
      call dsyev('V','U',n,ac,n,values,work,lwork,info)
      if(info/=0)return
      vectors=ac
      desc=.true.
      if(present(descending))desc=descending
      if(desc) then
         allocate(tmp(n))
         do i=1,n/2
            j=n-i+1
            tv=values(i)
            values(i)=values(j)
            values(j)=tv
            tmp=vectors(:,i)
            vectors(:,i)=vectors(:,j)
            vectors(:,j)=tmp
         end do
      end if
   end subroutine symmetric_eigen

   subroutine orthogonal_complement(a,b,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: b(:,:)
      integer, intent(out) :: info
      integer :: m,n,k,lwork,j
      real(dp), allocatable :: qmat(:,:),tau(:),work(:)
      real(dp) :: query(1)
      m=size(a,1)
      n=size(a,2)
      k=min(m,n)
      if(n>=m) then
      allocate(b(m,0))
      info=0
      return
      end if
      allocate(qmat(m,m),tau(k))
      qmat=0.0_dp
      qmat(:,1:n)=a
! Fill remaining columns before QR so storage is valid; dorgqr constructs Q.
      call dgeqrf(m,n,qmat,m,tau,query,-1,info)
      if(info/=0)return
      lwork=max(1,int(query(1)))
      allocate(work(lwork))
      qmat(:,1:n)=a
      call dgeqrf(m,n,qmat,m,tau,work,lwork,info)
      if(info/=0)return
      deallocate(work)
      call dorgqr(m,m,n,qmat,m,tau,query,-1,info)
      if(info/=0)return
      lwork=max(1,int(query(1)))
      allocate(work(lwork))
      call dorgqr(m,m,n,qmat,m,tau,work,lwork,info)
      if(info/=0)return
      allocate(b(m,m-n))
      b=qmat(:,n+1:m)
   end subroutine orthogonal_complement

   real(dp) function logdet_spd(a,info) result(v)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: l(:,:)
      integer :: i
      call chol_lower(a,l,info)
      if(info/=0) then
      v=-huge(1.0_dp)
      return
      end if
      v=0.0_dp
      do i=1,size(l,1)
      v=v+2.0_dp*log(l(i,i))
      end do
   end function logdet_spd

   real(dp) function determinant_spd(a,info) result(v)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out) :: info
      real(dp) :: ld
      ld=logdet_spd(a,info)
      if(info==0) then
      v=exp(ld)
      else
      v=0.0_dp
      end if
   end function determinant_spd
end module urca_linalg
