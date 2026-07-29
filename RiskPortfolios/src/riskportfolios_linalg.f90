! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_linalg
   use riskportfolios_kinds, only : dp
   implicit none
   private
   public :: solve_linear, symmetric_eigen, is_positive_definite

   interface
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         integer, intent(in) :: n, nrhs, lda, ldb
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
      end subroutine dgesv

      subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
         import :: dp
         character(len=1), intent(in) :: jobz, uplo
         integer, intent(in) :: n, lda, lwork
         integer, intent(out) :: info
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: w(*), work(*)
      end subroutine dsyev

      subroutine dpotrf(uplo, n, a, lda, info)
         import :: dp
         character(len=1), intent(in) :: uplo
         integer, intent(in) :: n, lda
         integer, intent(out) :: info
         real(dp), intent(inout) :: a(lda, *)
      end subroutine dpotrf
   end interface

contains

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: aa(:, :), bb(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n

      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n) then
         info = -1
         allocate(x(0))
         return
      end if
      allocate(aa(n, n), bb(n, 1), ipiv(n), x(n))
      aa = a
      bb(:, 1) = b
      call dgesv(n, 1, aa, n, ipiv, bb, n, info)
      if (info == 0) then
         x = bb(:, 1)
      else
         x = 0.0_dp
      end if
   end subroutine solve_linear

   subroutine symmetric_eigen(a, eigenvalues, eigenvectors, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: eigenvalues(:), eigenvectors(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: work(:)
      real(dp) :: work_query(1)
      integer :: n, lwork

      n = size(a, 1)
      if (size(a, 2) /= n) then
         info = -1
         allocate(eigenvalues(0), eigenvectors(0, 0))
         return
      end if
      allocate(eigenvalues(n), eigenvectors(n, n))
      eigenvectors = 0.5_dp * (a + transpose(a))
      lwork = -1
      call dsyev('V', 'U', n, eigenvectors, n, eigenvalues, work_query, lwork, info)
      if (info /= 0) return
      lwork = max(1, int(work_query(1)))
      allocate(work(lwork))
      call dsyev('V', 'U', n, eigenvectors, n, eigenvalues, work, lwork, info)
   end subroutine symmetric_eigen

   logical function is_positive_definite(a, tolerance) result(ok)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: aa(:, :)
      real(dp) :: tol
      integer :: n, info

      n = size(a, 1)
      if (size(a, 2) /= n) then
         ok = .false.
         return
      end if
      tol = 0.0_dp
      if (present(tolerance)) tol = tolerance
      allocate(aa(n, n))
      aa = 0.5_dp * (a + transpose(a))
      if (tol > 0.0_dp) aa = aa + tol * identity_matrix(n)
      call dpotrf('U', n, aa, n, info)
      ok = info == 0
   end function is_positive_definite

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

end module riskportfolios_linalg
