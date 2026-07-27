! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License version 2 only.
module pa_utils
  use pa_kinds, only: dp
  implicit none
  private
  public :: unique_int, average_ranks, quintile_groups, build_design_matrix
  public :: geometric_product

contains

  subroutine unique_int(x, values)
    integer, intent(in) :: x(:)
    integer, allocatable, intent(out) :: values(:)
    integer, allocatable :: tmp(:)
    integer :: i, n

    allocate(tmp(size(x)))
    n = 0
    do i = 1, size(x)
      if (n == 0 .or. .not. any(tmp(1:n) == x(i))) then
        n = n + 1
        tmp(n) = x(i)
      end if
    end do
    allocate(values(n))
    if (n > 0) values = tmp(1:n)
  end subroutine unique_int

  subroutine average_ranks(x, ranks)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: ranks(:)
    integer, allocatable :: idx(:)
    integer :: i, j, key, left, right
    real(dp) :: avg

    allocate(ranks(size(x)), idx(size(x)))
    do i = 1, size(x)
      idx(i) = i
    end do

    do i = 2, size(x)
      key = idx(i)
      j = i - 1
      do while (j >= 1)
        if (x(idx(j)) < x(key)) exit
        if (.not. (x(idx(j)) < x(key)) .and. .not. (x(key) < x(idx(j))) .and. idx(j) < key) exit
        idx(j + 1) = idx(j)
        j = j - 1
      end do
      idx(j + 1) = key
    end do

    left = 1
    do while (left <= size(x))
      right = left
      do while (right < size(x))
        if (x(idx(right + 1)) < x(idx(left)) .or. x(idx(left)) < x(idx(right + 1))) exit
        right = right + 1
      end do
      avg = 0.5_dp * real(left + right, dp)
      do i = left, right
        ranks(idx(i)) = avg
      end do
      left = right + 1
    end do
  end subroutine average_ranks

  subroutine quintile_groups(x, group)
    real(dp), intent(in) :: x(:)
    integer, allocatable, intent(out) :: group(:)
    real(dp), allocatable :: ranks(:)
    integer :: i, n

    n = size(x)
    allocate(group(n))
    if (n == 0) return
    call average_ranks(x, ranks)
    do i = 1, n
      group(i) = ceiling(5.0_dp * ranks(i) / real(n, dp))
      group(i) = max(1, min(5, group(i)))
    end do
  end subroutine quintile_groups

  subroutine build_design_matrix(categorical, n_levels, numeric, x, group_start, group_end, status)
    integer, intent(in) :: categorical(:, :)
    integer, intent(in) :: n_levels(:)
    real(dp), intent(in) :: numeric(:, :)
    real(dp), allocatable, intent(out) :: x(:, :)
    integer, allocatable, intent(out) :: group_start(:), group_end(:)
    integer, intent(out) :: status

    integer :: n, ncat, nnum, j, lev, col, first_level, ncol

    status = 0
    n = size(categorical, 1)
    ncat = size(categorical, 2)
    nnum = size(numeric, 2)
    if (size(numeric, 1) /= n .or. size(n_levels) /= ncat) then
      status = 1
      allocate(x(0, 0), group_start(0), group_end(0))
      return
    end if
    if (ncat > 0 .and. any(n_levels < 1)) then
      status = 2
      allocate(x(0, 0), group_start(0), group_end(0))
      return
    end if

    ncol = 0
    do j = 1, ncat
      if (j == 1) then
        ncol = ncol + n_levels(j)
      else
        ncol = ncol + max(0, n_levels(j) - 1)
      end if
    end do
    ncol = ncol + nnum

    allocate(x(n, ncol), group_start(ncat + nnum), group_end(ncat + nnum))
    x = 0.0_dp
    col = 0
    do j = 1, ncat
      group_start(j) = col + 1
      if (j == 1) then
        first_level = 1
      else
        first_level = 2
      end if
      do lev = first_level, n_levels(j)
        col = col + 1
        where (categorical(:, j) == lev) x(:, col) = 1.0_dp
      end do
      group_end(j) = col
      if (any(categorical(:, j) < 1) .or. any(categorical(:, j) > n_levels(j))) status = 3
    end do

    do j = 1, nnum
      col = col + 1
      x(:, col) = numeric(:, j)
      group_start(ncat + j) = col
      group_end(ncat + j) = col
    end do
  end subroutine build_design_matrix

  pure real(dp) function geometric_product(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 1.0_dp
    else
      value = product(1.0_dp + x)
    end if
  end function geometric_product

end module pa_utils
