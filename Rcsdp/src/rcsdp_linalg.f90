! BLAS/LAPACK interfaces used by the CSDP translation.
! See LICENSE (CPL-1.0).
module rcsdp_linalg
   use rcsdp_kinds, only : dp
   implicit none
   private
   public :: potrf_upper, potrs_upper, trtri_upper, symmetric_eigenvalues
   public :: solve_spd_factored, solve_spd, dot_product_dp, norm2_dp

   interface
      subroutine dpotrf(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dpotrf

      subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(in) :: a(lda,*)
         real(dp), intent(inout) :: b(ldb,*)
         integer, intent(out) :: info
      end subroutine dpotrs

      subroutine dtrtri(uplo, diag, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo, diag
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda,*)
         integer, intent(out) :: info
      end subroutine dtrtri

      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda,*)
         real(dp), intent(out) :: w(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dsyev
   end interface

contains

   subroutine potrf_upper(a, info)
      real(dp), intent(inout) :: a(:,:)
      integer, intent(out) :: info
      integer :: n, i, j
      n = size(a,1)
      if (size(a,2) /= n) then
         info = -1
         return
      end if
      call dpotrf('U', n, a, n, info)
      if (info == 0) then
         do j = 1, n - 1
            do i = j + 1, n
               a(i,j) = 0.0_dp
            end do
         end do
      end if
   end subroutine potrf_upper

   subroutine potrs_upper(a, b, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(inout) :: b(:)
      integer, intent(out) :: info
      real(dp), allocatable :: rhs(:,:)
      integer :: n
      n = size(a,1)
      if (size(a,2) /= n .or. size(b) /= n) then
         info = -1
         return
      end if
      allocate(rhs(n,1))
      rhs(:,1) = b
      call dpotrs('U', n, 1, a, n, rhs, n, info)
      if (info == 0) b = rhs(:,1)
   end subroutine potrs_upper

   subroutine trtri_upper(a, info)
      real(dp), intent(inout) :: a(:,:)
      integer, intent(out) :: info
      integer :: n
      n = size(a,1)
      if (size(a,2) /= n) then
         info = -1
         return
      end if
      call dtrtri('U', 'N', n, a, n, info)
   end subroutine trtri_upper

   subroutine symmetric_eigenvalues(a, w, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), allocatable :: acopy(:,:), work(:)
      real(dp) :: workq(1)
      integer :: n, lwork
      n = size(a,1)
      if (size(a,2) /= n) then
         info = -1
         allocate(w(0))
         return
      end if
      allocate(acopy(n,n), w(n))
      acopy = 0.5_dp*(a + transpose(a))
      lwork = -1
      call dsyev('N', 'U', n, acopy, n, w, workq, lwork, info)
      if (info /= 0) return
      lwork = max(1, int(workq(1)))
      allocate(work(lwork))
      call dsyev('N', 'U', n, acopy, n, w, work, lwork, info)
   end subroutine symmetric_eigenvalues

   subroutine solve_spd_factored(r, b, info)
      real(dp), intent(in) :: r(:,:)
      real(dp), intent(inout) :: b(:)
      integer, intent(out) :: info
      call potrs_upper(r, b, info)
   end subroutine solve_spd_factored

   subroutine solve_spd(a, b, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(inout) :: b(:)
      integer, intent(out) :: info
      real(dp), allocatable :: r(:,:)
      allocate(r(size(a,1),size(a,2)))
      r = a
      call potrf_upper(r, info)
      if (info == 0) call solve_spd_factored(r, b, info)
   end subroutine solve_spd

   pure real(dp) function dot_product_dp(x,y) result(v)
      real(dp), intent(in) :: x(:), y(:)
      v = dot_product(x,y)
   end function dot_product_dp

   pure real(dp) function norm2_dp(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(max(0.0_dp, dot_product(x,x)))
   end function norm2_dp

end module rcsdp_linalg
