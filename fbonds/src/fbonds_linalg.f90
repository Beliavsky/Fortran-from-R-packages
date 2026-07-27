! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
module fbonds_linalg
   use fbonds_kinds, only : dp
   implicit none
   private
   public :: least_squares

   interface
      subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
         import dp
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         integer, intent(out) :: rank, info
         real(dp), intent(inout) :: a(lda, *), b(ldb, *)
         real(dp), intent(out) :: s(*)
         real(dp), intent(in) :: rcond
         real(dp), intent(inout) :: work(*)
      end subroutine dgelss
   end interface
contains
   subroutine least_squares(x, y, beta, rank, info)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), intent(out) :: beta(:)
      integer, intent(out), optional :: rank, info
      real(dp), allocatable :: a(:, :), b(:, :), s(:), work(:)
      real(dp) :: work_query(1)
      integer :: m, n, lda, ldb, irank, iinfo, lwork

      m = size(x, 1)
      n = size(x, 2)
      if (size(y) /= m .or. size(beta) /= n .or. m < 1 .or. n < 1) then
         beta = 0.0_dp
         if (present(rank)) rank = 0
         if (present(info)) info = -1
         return
      end if

      lda = max(1, m)
      ldb = max(m, n)
      allocate(a(lda, n), b(ldb, 1), s(min(m, n)))
      a = x
      b = 0.0_dp
      b(1:m, 1) = y
      lwork = -1
      call dgelss(m, n, 1, a, lda, b, ldb, s, -1.0_dp, irank, work_query, lwork, iinfo)
      if (iinfo /= 0) then
         beta = 0.0_dp
         if (present(rank)) rank = 0
         if (present(info)) info = iinfo
         return
      end if
      lwork = max(1, ceiling(work_query(1)))
      allocate(work(lwork))
      a = x
      b = 0.0_dp
      b(1:m, 1) = y
      call dgelss(m, n, 1, a, lda, b, ldb, s, -1.0_dp, irank, work, lwork, iinfo)
      if (iinfo == 0) then
         beta = b(1:n, 1)
      else
         beta = 0.0_dp
      end if
      if (present(rank)) rank = irank
      if (present(info)) info = iinfo
   end subroutine least_squares
end module fbonds_linalg
