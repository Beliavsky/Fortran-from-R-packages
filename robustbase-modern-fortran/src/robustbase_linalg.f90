! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_linalg
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: solve_linear, least_squares, symmetric_eigen, invert_symmetric, matrix_rank, covariance_matrix
   interface
      subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
         import dp
         integer, intent(in) :: n,nrhs,lda,ldb
         integer, intent(out) :: ipiv(*)
         real(dp), intent(inout) :: a(lda,*),b(ldb,*)
         integer, intent(out) :: info
      end subroutine dgesv
      subroutine dgels(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)
         import dp
         character(len=1), intent(in) :: trans
         integer, intent(in) :: m,n,nrhs,lda,ldb,lwork
         real(dp), intent(inout) :: a(lda,*),b(ldb,*),work(*)
         integer, intent(out) :: info
      end subroutine dgels
      subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         import dp
         character(len=1), intent(in) :: jobz,uplo
         integer, intent(in) :: n,lda,lwork
         real(dp), intent(inout) :: a(lda,*),work(*)
         real(dp), intent(out) :: w(*)
         integer, intent(out) :: info
      end subroutine dsyev
   end interface
contains
   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:), bb(:,:)
      integer, allocatable :: ipiv(:)
      integer :: n
      n = size(a,1)
      if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) error stop "solve_linear: size mismatch"
      aa = a
      allocate(bb(n,1), ipiv(n))
      bb(:,1) = b
      call dgesv(n,1,aa,n,ipiv,bb,n,info)
      if (info == 0) x = bb(:,1)
   end subroutine solve_linear

   subroutine least_squares(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:), bb(:,:), work(:)
      real(dp) :: wq(1)
      integer :: m,n,lwork,ldb
      m=size(a,1); n=size(a,2); ldb=max(m,n)
      if(size(b)/=m .or. size(x)/=n) error stop "least_squares: size mismatch"
      aa=a
      allocate(bb(ldb,1)); bb=0.0_dp; bb(1:m,1)=b
      call dgels('N',m,n,1,aa,m,bb,ldb,wq,-1,info)
      if(info/=0) return
      lwork=max(1,int(wq(1)))
      allocate(work(lwork)); aa=a; bb=0.0_dp; bb(1:m,1)=b
      call dgels('N',m,n,1,aa,m,bb,ldb,work,lwork,info)
      if(info==0) x=bb(1:n,1)
   end subroutine least_squares

   subroutine symmetric_eigen(a, values, vectors, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: values(:), vectors(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:), work(:)
      real(dp) :: wq(1)
      integer :: n,lwork
      n=size(a,1)
      if(size(a,2)/=n .or. size(values)/=n .or. any(shape(vectors)/=[n,n])) error stop "symmetric_eigen: size mismatch"
      aa=a
      call dsyev('V','U',n,aa,n,values,wq,-1,info)
      if(info/=0) return
      lwork=max(1,int(wq(1)))
      allocate(work(lwork)); aa=a
      call dsyev('V','U',n,aa,n,values,work,lwork,info)
      if(info==0) vectors=aa
   end subroutine symmetric_eigen

   subroutine invert_symmetric(a, ainv, info, ridge)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: vals(:), vecs(:,:)
      real(dp) :: tol, rr
      integer :: n,i
      n=size(a,1); rr=0.0_dp; if(present(ridge)) rr=ridge
      allocate(vals(n),vecs(n,n))
      call symmetric_eigen(a,vals,vecs,info); if(info/=0) return
      tol=max(1.0e-12_dp,epsilon(1.0_dp)*maxval(abs(vals))*real(n,dp))
      ainv=0.0_dp
      do i=1,n
         if(vals(i)+rr > tol) ainv=ainv+outer(vecs(:,i),vecs(:,i))/(vals(i)+rr)
      end do
   contains
      pure function outer(x,y) result(z)
         real(dp),intent(in)::x(:),y(:)
         real(dp)::z(size(x),size(y))
         integer::j
         do j=1,size(y); z(:,j)=x*y(j); end do
      end function outer
   end subroutine invert_symmetric

   integer function matrix_rank(a, tol) result(r)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: gram(:,:), vals(:), vecs(:,:)
      real(dp) :: t
      integer :: info,n
      gram=matmul(transpose(a),a); n=size(gram,1)
      allocate(vals(n),vecs(n,n)); call symmetric_eigen(gram,vals,vecs,info)
      if(info/=0) then; r=0; return; end if
      if(present(tol)) then; t=tol; else; t=max(1.0e-12_dp,epsilon(1.0_dp)*maxval(abs(vals))*real(max(size(a,1),size(a,2)),dp)); end if
      r=count(vals>t*t)
   end function matrix_rank

   subroutine covariance_matrix(x, center, cov)
      real(dp), intent(in) :: x(:,:), center(:)
      real(dp), intent(out) :: cov(:,:)
      real(dp), allocatable :: z(:,:)
      integer :: n,j
      n=size(x,1)
      if(size(center)/=size(x,2) .or. any(shape(cov)/=[size(x,2),size(x,2)])) error stop "covariance_matrix: size mismatch"
      allocate(z(n,size(x,2)))
      do j=1,size(x,2); z(:,j)=x(:,j)-center(j); end do
      if(n>1) then; cov=matmul(transpose(z),z)/real(n-1,dp); else; cov=0.0_dp; end if
   end subroutine covariance_matrix
end module robustbase_linalg
