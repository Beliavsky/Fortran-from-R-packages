! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_preprocess
  use fextremes_kinds, only: dp
  use fextremes_stats, only: sort_real
  implicit none
  private
  public :: point_process_result, decluster_result
  public :: block_maxima, find_threshold, point_process, decluster

  type :: point_process_result
    real(dp), allocatable :: values(:)
    integer, allocatable :: indices(:)
  end type point_process_result

  type :: decluster_result
    real(dp), allocatable :: maxima(:)
    integer, allocatable :: max_indices(:)
    integer, allocatable :: from_indices(:)
    integer, allocatable :: to_indices(:)
  end type decluster_result
contains
  subroutine block_maxima(x, block, maxima, max_indices)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: block
    real(dp), allocatable, intent(out) :: maxima(:)
    integer, allocatable, intent(out) :: max_indices(:)
    integer :: nb, b, lo, hi, j, loc
    if (block <= 0 .or. size(x) == 0) then
      allocate(maxima(0), max_indices(0)); return
    end if
    nb = (size(x) + block - 1) / block
    allocate(maxima(nb), max_indices(nb))
    do b = 1, nb
      lo = (b - 1) * block + 1
      hi = min(size(x), b * block)
      loc = lo
      do j = lo + 1, hi
        if (x(j) > x(loc)) loc = j
      end do
      maxima(b) = x(loc)
      max_indices(b) = loc
    end do
  end subroutine block_maxima

  real(dp) function find_threshold(x, n_extremes) result(u)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n_extremes
    real(dp), allocatable :: work(:), unique_vals(:)
    integer :: n, k, m, i
    n = size(x)
    if (n == 0) then
      u = 0.0_dp; return
    end if
    allocate(work(n), unique_vals(n)); work = x
    call sort_real(work, ascending=.false.)
    k = max(1, min(n, n_extremes))
    m = 1; unique_vals(1) = work(1)
    do i = 2, n
      if (abs(work(i)-unique_vals(m)) > 0.0_dp) then
        m = m + 1; unique_vals(m) = work(i)
      end if
    end do
    i = 1
    do while (i <= m)
      if (abs(unique_vals(i)-work(k)) <= 0.0_dp) exit
      i = i + 1
    end do
    u = unique_vals(min(m, i + 1))
  end function find_threshold

  subroutine point_process(x, threshold, result)
    real(dp), intent(in) :: x(:), threshold
    type(point_process_result), intent(out) :: result
    integer :: n, i, j
    n = count(x > threshold)
    allocate(result%values(n), result%indices(n))
    j = 0
    do i = 1, size(x)
      if (x(i) > threshold) then
        j = j + 1; result%values(j) = x(i); result%indices(j) = i
      end if
    end do
  end subroutine point_process

  subroutine decluster(values, positions, run_length, result)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: positions(:), run_length
    type(decluster_result), intent(out) :: result
    integer :: ncl, i, j, start_i, end_i, max_i
    if (size(values) /= size(positions) .or. size(values) == 0) then
      allocate(result%maxima(0), result%max_indices(0), result%from_indices(0), result%to_indices(0)); return
    end if
    ncl = 1
    do i = 2, size(values)
      if (positions(i) - positions(i - 1) > run_length) ncl = ncl + 1
    end do
    allocate(result%maxima(ncl), result%max_indices(ncl), result%from_indices(ncl), result%to_indices(ncl))
    start_i = 1; j = 0
    do i = 2, size(values)
      if (positions(i) - positions(i - 1) > run_length) then
        end_i = i - 1
        j = j + 1
        max_i = start_i
        if (end_i > start_i) max_i = start_i - 1 + &
          maxloc(values(start_i:end_i), dim=1)
        result%maxima(j) = values(max_i)
        result%max_indices(j) = positions(max_i)
        result%from_indices(j) = positions(start_i)
        result%to_indices(j) = positions(end_i)
        start_i = i
      end if
    end do
    end_i = size(values)
    j = j + 1
    max_i = start_i
    if (end_i > start_i) max_i = start_i - 1 + &
      maxloc(values(start_i:end_i), dim=1)
    result%maxima(j) = values(max_i)
    result%max_indices(j) = positions(max_i)
    result%from_indices(j) = positions(start_i)
    result%to_indices(j) = positions(end_i)
  end subroutine decluster
end module fextremes_preprocess
