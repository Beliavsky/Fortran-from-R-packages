! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_linalg
   use garchx_kinds, only : dp
   implicit none
   private
   public :: invert_matrix, solve_linear, cholesky_lower, matrix_rank

   interface
      subroutine dgetrf(m, n, a, lda, ipiv, info)
         import dp
         integer, intent(in) :: m, n, lda
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine dgetrf
      subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
         import dp
         integer, intent(in) :: n, lda, lwork
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         real(dp), intent(inout) :: work(*)
         integer, intent(out) :: info
      end subroutine dgetri
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import dp
         integer, intent(in) :: n, nrhs, lda, ldb
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
      end subroutine dgesv
      subroutine dpotrf(uplo, n, a, lda, info)
         import dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(out) :: info
      end subroutine dpotrf
      subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
         import dp
         character(len=1), intent(in) :: jobu, jobvt
         integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
         integer, intent(out) :: info
      end subroutine dgesvd
   end interface
contains
   subroutine invert_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: status
      integer :: n, info, lwork
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: work(:)

      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         status = 1
         allocate(ainv(0, 0))
         return
      end if
      allocate(ainv(n, n), ipiv(n))
      ainv = a
      call dgetrf(n, n, ainv, n, ipiv, info)
      if (info /= 0) then
         status = info
         ainv = 0.0_dp
         return
      end if
      lwork = max(1, 64*n)
      allocate(work(lwork))
      call dgetri(n, ainv, n, ipiv, work, lwork, info)
      status = info
      if (info /= 0) ainv = 0.0_dp
   end subroutine invert_matrix

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: status
      integer :: n, info
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: aa(:, :), bb(:, :)

      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. n < 1) then
         status = 1
         allocate(x(0))
         return
      end if
      allocate(aa(n, n), bb(n, 1), ipiv(n), x(n))
      aa = a
      bb(:, 1) = b
      call dgesv(n, 1, aa, n, ipiv, bb, n, info)
      status = info
      if (info == 0) then
         x = bb(:, 1)
      else
         x = 0.0_dp
      end if
   end subroutine solve_linear

   subroutine cholesky_lower(a, l, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: l(:, :)
      integer, intent(out) :: status
      integer :: n, info, i

      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         status = 1
         allocate(l(0, 0))
         return
      end if
      allocate(l(n, n))
      l = a
      call dpotrf('L', n, l, n, info)
      status = info
      if (info /= 0) then
         l = 0.0_dp
         return
      end if
      do i = 1, n
         if (i < n) l(i, i+1:n) = 0.0_dp
      end do
   end subroutine cholesky_lower

   integer function matrix_rank(a, tol) result(rank_value)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tol
      integer :: m, n, k, info, lwork
      real(dp) :: threshold
      real(dp), allocatable :: aa(:, :), s(:), u(:, :), vt(:, :), work(:)

      m = size(a, 1)
      n = size(a, 2)
      k = min(m, n)
      if (m == 0 .or. n == 0) then
         rank_value = 0
         return
      end if
      allocate(aa(m, n), s(k), u(1, 1), vt(1, 1))
      aa = a
      lwork = max(1, 8*max(m, n))
      allocate(work(lwork))
      call dgesvd('N', 'N', m, n, aa, m, s, u, 1, vt, 1, work, lwork, info)
      if (info /= 0) then
         rank_value = 0
         return
      end if
      if (present(tol)) then
         threshold = tol
      else
         threshold = max(m, n)*epsilon(1.0_dp)*max(1.0_dp, s(1))
      end if
      rank_value = count(s > threshold)
   end function matrix_rank
end module garchx_linalg
