! SPDX-License-Identifier: GPL-2.0-or-later
! Thin ABI compatibility wrappers for the legacy mclust BLAS/LAPACK calls.
subroutine daxpy(n,a,x,incx,y,incy)
  use la_blas_d, only: backend => la_daxpy
  integer :: n,incx,incy
  double precision :: a,x(*),y(*)
  call backend(n,a,x,incx,y,incy)
end subroutine daxpy

subroutine dcopy(n,x,incx,y,incy)
  use la_blas_d, only: backend => la_dcopy
  integer :: n,incx,incy
  double precision :: x(*),y(*)
  call backend(n,x,incx,y,incy)
end subroutine dcopy

double precision function ddot(n,x,incx,y,incy)
  use la_blas_d, only: backend => la_ddot
  integer :: n,incx,incy
  double precision :: x(*),y(*)
  ddot=backend(n,x,incx,y,incy)
end function ddot

subroutine dgemm(ta,tb,m,n,k,a,x,ldx,y,ldy,b,z,ldz)
  use la_blas_d, only: backend => la_dgemm
  character(len=1) :: ta,tb
  integer :: m,n,k,ldx,ldy,ldz
  double precision :: a,b,x(ldx,*),y(ldy,*),z(ldz,*)
  call backend(ta,tb,m,n,k,a,x,ldx,y,ldy,b,z,ldz)
end subroutine dgemm

subroutine dgemv(t,m,n,a,x,ldx,y,incy,b,z,incz)
  use la_blas_d, only: backend => la_dgemv
  character(len=1) :: t
  integer :: m,n,ldx,incy,incz
  double precision :: a,b,x(ldx,*),y(*),z(*)
  call backend(t,m,n,a,x,ldx,y,incy,b,z,incz)
end subroutine dgemv

subroutine dger(m,n,a,x,incx,y,incy,z,ldz)
  use la_blas_d, only: backend => la_dger
  integer :: m,n,incx,incy,ldz
  double precision :: a,x(*),y(*),z(ldz,*)
  call backend(m,n,a,x,incx,y,incy,z,ldz)
end subroutine dger

subroutine drot(n,x,incx,y,incy,c,s)
  use la_blas_d, only: backend => la_drot
  integer :: n,incx,incy
  double precision :: x(*),y(*),c,s
  call backend(n,x,incx,y,incy,c,s)
end subroutine drot

subroutine drotg(a,b,c,s)
  use la_blas_d, only: backend => la_drotg
  double precision :: a,b,c,s
  call backend(a,b,c,s)
end subroutine drotg

subroutine dscal(n,a,x,incx)
  use la_blas_d, only: backend => la_dscal
  integer :: n,incx
  double precision :: a,x(*)
  call backend(n,a,x,incx)
end subroutine dscal

subroutine dswap(n,x,incx,y,incy)
  use la_blas_d, only: backend => la_dswap
  integer :: n,incx,incy
  double precision :: x(*),y(*)
  call backend(n,x,incx,y,incy)
end subroutine dswap

subroutine dsyrk(u,t,n,k,a,x,ldx,b,y,ldy)
  use la_blas_d, only: backend => la_dsyrk
  character(len=1) :: u,t
  integer :: n,k,ldx,ldy
  double precision :: a,b,x(ldx,*),y(ldy,*)
  call backend(u,t,n,k,a,x,ldx,b,y,ldy)
end subroutine dsyrk

subroutine dtrsv(u,t,d,n,a,lda,x,incx)
  use la_blas_d, only: backend => la_dtrsv
  character(len=1) :: u,t,d
  integer :: n,lda,incx
  double precision :: a(lda,*),x(*)
  call backend(u,t,d,n,a,lda,x,incx)
end subroutine dtrsv

subroutine dgesvd(ju,jv,m,n,a,lda,s,u,ldu,vt,ldvt,w,lw,info)
  use la_lapack_d, only: backend => la_dgesvd
  character(len=1) :: ju,jv
  integer :: m,n,lda,ldu,ldvt,lw,info
  double precision :: a(lda,*),s(*),u(ldu,*),vt(ldvt,*),w(*)
  call backend(ju,jv,m,n,a,lda,s,u,ldu,vt,ldvt,w,lw,info)
end subroutine dgesvd

subroutine dlassq(n,x,incx,scale,sumsq)
  use la_lapack_d, only: backend => la_dlassq
  integer :: n,incx
  double precision :: x(*),scale,sumsq
  call backend(n,x,incx,scale,sumsq)
end subroutine dlassq

subroutine dpotrf(u,n,a,lda,info)
  use la_lapack_d, only: backend => la_dpotrf
  character(len=1) :: u
  integer :: n,lda,info
  double precision :: a(lda,*)
  call backend(u,n,a,lda,info)
end subroutine dpotrf

subroutine dpotri(u,n,a,lda,info)
  use la_lapack_d, only: backend => la_dpotri
  character(len=1) :: u
  integer :: n,lda,info
  double precision :: a(lda,*)
  call backend(u,n,a,lda,info)
end subroutine dpotri

subroutine dsyev(j,u,n,a,lda,w,work,lwork,info)
  use la_lapack_d, only: backend => la_dsyev
  character(len=1) :: j,u
  integer :: n,lda,lwork,info
  double precision :: a(lda,*),w(*),work(*)
  call backend(j,u,n,a,lda,w,work,lwork,info)
end subroutine dsyev

subroutine dsyevd(j,u,n,a,lda,w,work,lwork,iwork,liwork,info)
  use la_lapack_d, only: backend => la_dsyevd
  character(len=1) :: j,u
  integer :: n,lda,lwork,liwork,info,iwork(*)
  double precision :: a(lda,*),w(*),work(*)
  call backend(j,u,n,a,lda,w,work,lwork,iwork,liwork,info)
end subroutine dsyevd

subroutine dsyevx(j,r,u,n,a,lda,vl,vu,il,iu,tol,m,w,z,ldz,work,lwork,iwork,ifail,info)
  use la_lapack_d, only: backend => la_dsyevx
  character(len=1) :: j,r,u
  integer :: n,lda,il,iu,m,ldz,lwork,info,iwork(*),ifail(*)
  double precision :: a(lda,*),vl,vu,tol,w(*),z(ldz,*),work(*)
  call backend(j,r,u,n,a,lda,vl,vu,il,iu,tol,m,w,z,ldz,work,lwork,iwork,ifail,info)
end subroutine dsyevx
