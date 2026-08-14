! SPDX-License-Identifier: GPL-2.0-only
module ks_lapack
   use ks_kinds, only: dp
   implicit none
   private
   public :: dpotrf, dpotri, dpotrs, dsyev, dgesv
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
      subroutine dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(in) :: a(lda,*)
         real(dp), intent(inout) :: b(ldb,*)
         integer, intent(out) :: info
      end subroutine dpotrs
      subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsyev
      subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
         integer, intent(out) :: ipiv(*), info
      end subroutine dgesv
   end interface
end module ks_lapack
