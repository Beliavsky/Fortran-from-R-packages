! SPDX-License-Identifier: GPL-2.0-or-later
module nleqslv_blas_lapack
   implicit none
   private
   public :: daxpy, dcopy, ddot, dgemv, dgeqrf, dlacpy, dlamch, dlantr
   public :: dlartg, dnrm2, dorgqr, dormqr, drot, dscal, dtrcon, dtrmm
   public :: dtrmv, dtrsv, idamax

   interface
      subroutine daxpy(n, da, dx, incx, dy, incy)
         integer :: n, incx, incy
         double precision :: da, dx(*), dy(*)
      end subroutine daxpy

      subroutine dcopy(n, dx, incx, dy, incy)
         integer :: n, incx, incy
         double precision :: dx(*), dy(*)
      end subroutine dcopy

      double precision function ddot(n, dx, incx, dy, incy)
         integer :: n, incx, incy
         double precision :: dx(*), dy(*)
      end function ddot

      subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         character(len=1) :: trans
         integer :: m, n, lda, incx, incy
         double precision :: alpha, beta, a(lda,*), x(*), y(*)
      end subroutine dgemv

      subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
         integer :: m, n, lda, lwork, info
         double precision :: a(lda,*), tau(*), work(*)
      end subroutine dgeqrf

      subroutine dlacpy(uplo, m, n, a, lda, b, ldb)
         character(len=1) :: uplo
         integer :: m, n, lda, ldb
         double precision :: a(lda,*), b(ldb,*)
      end subroutine dlacpy

      double precision function dlamch(cmach)
         character(len=1) :: cmach
      end function dlamch

      double precision function dlantr(norm, uplo, diag, m, n, a, lda, work)
         character(len=1) :: norm, uplo, diag
         integer :: m, n, lda
         double precision :: a(lda,*), work(*)
      end function dlantr

      subroutine dlartg(f, g, cs, sn, r)
         double precision :: f, g, cs, sn, r
      end subroutine dlartg

      double precision function dnrm2(n, x, incx)
         integer :: n, incx
         double precision :: x(*)
      end function dnrm2

      subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
         integer :: m, n, k, lda, lwork, info
         double precision :: a(lda,*), tau(*), work(*)
      end subroutine dorgqr

      subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
         character(len=1) :: side, trans
         integer :: m, n, k, lda, ldc, lwork, info
         double precision :: a(lda,*), tau(*), c(ldc,*), work(*)
      end subroutine dormqr

      subroutine drot(n, dx, incx, dy, incy, c, s)
         integer :: n, incx, incy
         double precision :: dx(*), dy(*), c, s
      end subroutine drot

      subroutine dscal(n, da, dx, incx)
         integer :: n, incx
         double precision :: da, dx(*)
      end subroutine dscal

      subroutine dtrcon(norm, uplo, diag, n, a, lda, rcond, work, iwork, info)
         character(len=1) :: norm, uplo, diag
         integer :: n, lda, iwork(*), info
         double precision :: a(lda,*), rcond, work(*)
      end subroutine dtrcon

      subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
         character(len=1) :: side, uplo, transa, diag
         integer :: m, n, lda, ldb
         double precision :: alpha, a(lda,*), b(ldb,*)
      end subroutine dtrmm

      subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
         character(len=1) :: uplo, trans, diag
         integer :: n, lda, incx
         double precision :: a(lda,*), x(*)
      end subroutine dtrmv

      subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
         character(len=1) :: uplo, trans, diag
         integer :: n, lda, incx
         double precision :: a(lda,*), x(*)
      end subroutine dtrsv

      integer function idamax(n, dx, incx)
         integer :: n, incx
         double precision :: dx(*)
      end function idamax
   end interface
end module nleqslv_blas_lapack
