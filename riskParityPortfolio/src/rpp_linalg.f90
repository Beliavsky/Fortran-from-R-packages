! SPDX-License-Identifier: GPL-3.0-only
! Linear-algebra helpers for the riskParityPortfolio Fortran port.
module rpp_linalg
   use rpp_kinds, only: dp
   implicit none
   private
   public :: solve_linear, solve_spd, vector_norm2, max_abs, matrix_trace
   public :: identity_matrix, symmetric_part

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*)
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
         integer, intent(out) :: info
      end subroutine dgesv
      subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
         integer, intent(out) :: info
      end subroutine dposv
   end interface
contains
   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: ac(:, :), bc(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(ac(n, n), bc(n, 1), ipiv(n))
      ac = a
      bc(:, 1) = b
      call dgesv(n, 1, ac, n, ipiv, bc, n, info)
      if (info == 0) x = bc(:, 1)
   end subroutine solve_linear

   subroutine solve_spd(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: ac(:, :), bc(:, :)
      integer :: n
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      allocate(ac(n, n), bc(n, 1))
      ac = a
      bc(:, 1) = b
      call dposv('U', n, 1, ac, n, bc, n, info)
      if (info == 0) x = bc(:, 1)
   end subroutine solve_spd

   pure real(dp) function vector_norm2(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(max(0.0_dp, dot_product(x, x)))
   end function vector_norm2

   pure real(dp) function max_abs(x) result(v)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         v = 0.0_dp
      else
         v = maxval(abs(x))
      end if
   end function max_abs

   pure real(dp) function matrix_trace(a) result(v)
      real(dp), intent(in) :: a(:, :)
      integer :: i
      v = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         v = v + a(i, i)
      end do
   end function matrix_trace

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure function symmetric_part(a) result(s)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: s(size(a, 1), size(a, 2))
      s = 0.5_dp * (a + transpose(a))
   end function symmetric_part
end module rpp_linalg
