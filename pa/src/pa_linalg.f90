! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 only.
module pa_linalg
  use pa_kinds, only: dp
  implicit none
  private
  public :: least_squares

  interface
    subroutine dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
      import :: dp
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(inout) :: b(ldb, *)
      integer, intent(inout) :: jpvt(*)
      real(dp), intent(in) :: rcond
      integer, intent(out) :: rank, info
      real(dp), intent(inout) :: work(*)
    end subroutine dgelsy
  end interface

contains

  subroutine least_squares(x, y, beta, rank, status)
    real(dp), intent(in) :: x(:, :)
    real(dp), intent(in) :: y(:)
    real(dp), allocatable, intent(out) :: beta(:)
    integer, intent(out) :: rank
    integer, intent(out) :: status

    real(dp), allocatable :: a(:, :), b(:, :), work(:)
    integer, allocatable :: jpvt(:)
    real(dp) :: work_query(1), rcond
    integer :: m, n, ldb, lwork, info

    m = size(x, 1)
    n = size(x, 2)
    rank = 0
    status = 0
    allocate(beta(n))
    beta = 0.0_dp

    if (size(y) /= m .or. m < 1 .or. n < 1) then
      status = 1
      return
    end if

    ldb = max(m, n)
    allocate(a(m, n), b(ldb, 1), jpvt(n))
    a = x
    b = 0.0_dp
    b(1:m, 1) = y
    jpvt = 0
    rcond = sqrt(epsilon(1.0_dp))

    lwork = -1
    call dgelsy(m, n, 1, a, m, b, ldb, jpvt, rcond, rank, work_query, lwork, info)
    if (info /= 0) then
      status = 2
      return
    end if

    lwork = max(1, ceiling(work_query(1)))
    allocate(work(lwork))
    a = x
    b = 0.0_dp
    b(1:m, 1) = y
    jpvt = 0
    call dgelsy(m, n, 1, a, m, b, ldb, jpvt, rcond, rank, work, lwork, info)
    if (info /= 0) then
      status = 3
      return
    end if
    beta = b(1:n, 1)
  end subroutine least_squares

end module pa_linalg
