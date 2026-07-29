! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_linalg
   use cccp_kinds, only : dp
   implicit none
   private
   public :: solve_system, symmetric_eigenvalues, spd_inverse_logdet
   public :: equality_particular, vector_norm2, symmetrize

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*)
         real(dp), intent(inout) :: a(lda,*), b(ldb,*)
         integer, intent(out) :: info
      end subroutine dgesv
      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*), work(*)
         integer, intent(out) :: info
      end subroutine dsyev
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
   end interface

contains

   pure function vector_norm2(x) result(ans)
      real(dp), intent(in) :: x(:)
      real(dp) :: ans
      ans = sqrt(max(0.0_dp, dot_product(x, x)))
   end function vector_norm2

   pure function symmetrize(a) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: s(size(a,1), size(a,2))
      s = 0.5_dp * (a + transpose(a))
   end function symmetrize

   subroutine solve_system(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:), bb(:,:)
      integer, allocatable :: ipiv(:)
      integer :: n, i

      n = size(b)
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(aa(n,n), bb(n,1), ipiv(n))
      aa = a
      bb(:,1) = b
      do i = 1, n
         aa(i,i) = aa(i,i) + 1.0e-12_dp * max(1.0_dp, abs(aa(i,i)))
      end do
      call dgesv(n, 1, aa, n, ipiv, bb, n, info)
      if (info == 0) x = bb(:,1)
   end subroutine solve_system

   subroutine symmetric_eigenvalues(a, w, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:,:), work(:)
      integer :: n, lwork

      n = size(a,1)
      if (size(a,2) /= n .or. size(w) /= n) then
         info = -1
         return
      end if
      lwork = max(1, 3*n - 1)
      allocate(aa(n,n), work(lwork))
      aa = symmetrize(a)
      call dsyev('N', 'U', n, aa, n, w, work, lwork, info)
   end subroutine symmetric_eigenvalues

   subroutine spd_inverse_logdet(a, ainv, logdet, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      real(dp), intent(out) :: logdet
      integer, intent(out) :: info
      integer :: n, i, j

      n = size(a,1)
      if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
         info = -1
         logdet = huge(1.0_dp)
         return
      end if
      ainv = symmetrize(a)
      call dpotrf('L', n, ainv, n, info)
      if (info /= 0) then
         logdet = huge(1.0_dp)
         return
      end if
      logdet = 0.0_dp
      do i = 1, n
         logdet = logdet + 2.0_dp * log(ainv(i,i))
      end do
      call dpotri('L', n, ainv, n, info)
      if (info /= 0) return
      do j = 1, n
         do i = 1, j - 1
            ainv(i,j) = ainv(j,i)
         end do
      end do
   end subroutine spd_inverse_logdet

   subroutine equality_particular(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aat(:,:), y(:)
      integer :: p, n

      p = size(a,1)
      n = size(a,2)
      if (size(b) /= p .or. size(x) /= n) then
         info = -1
         return
      end if
      if (p == 0) then
         x = 0.0_dp
         info = 0
         return
      end if
      allocate(aat(p,p), y(p))
      aat = matmul(a, transpose(a))
      call solve_system(aat, b, y, info)
      if (info == 0) x = matmul(transpose(a), y)
   end subroutine equality_particular

end module cccp_linalg
