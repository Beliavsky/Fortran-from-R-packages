! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
program analyze_csv
  use tserieschaos, only : dp, average_mutual_information, correlation_integral, &
    false_nearest_curve, lyapunov_stretching, lyapunov_linear_fit
  implicit none
  character(len=512) :: filename, line
  character(len=16) :: search_method, method_used
  real(dp), allocatable :: series(:), ami(:), fractions(:), stretching(:)
  integer, allocatable :: totals(:)
  real(dp) :: c2, intercept, exponent, value
  integer :: unit, ios, n, capacity, status, refs_used, comma
  logical :: first_line

  if (command_argument_count() < 1) then
    write(*, '(a)') "usage: analyze_csv FILE [auto|box|direct]"
    write(*, '(a)') "FILE may contain one numeric column or Date,Value rows with an optional header."
    stop 2
  end if
  call get_command_argument(1, filename)
  search_method = "auto"
  if (command_argument_count() >= 2) call get_command_argument(2, search_method)
  capacity = 1024
  allocate(series(capacity))
  n = 0
  first_line = .true.
  open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) error stop "cannot open CSV file"
  do
    read(unit, '(a)', iostat=ios) line
    if (ios < 0) exit
    if (ios > 0) error stop "error reading CSV file"
    comma = index(line, ',')
    if (comma > 0) then
      read(line(comma + 1:), *, iostat=ios) value
    else
      read(line, *, iostat=ios) value
    end if
    if (ios /= 0) then
      if (first_line) then
        first_line = .false.
        cycle
      end if
      error stop "non-numeric CSV value"
    end if
    first_line = .false.
    n = n + 1
    if (n > capacity) call grow(series, capacity)
    series(n) = value
  end do
  close(unit)
  if (n < 50) error stop "at least 50 observations are required"
  series = series(:n)

  call average_mutual_information(series, 16, 20, ami, status)
  call require(status == 0, "AMI calculation failed")
  call correlation_integral(series, 3, 1, 5, 0.1_dp * sample_sd(series), c2, status)
  call require(status == 0, "correlation integral failed")
  call false_nearest_curve(series, 5, 1, 5, 10.0_dp, 0.2_dp * sample_sd(series), &
    fractions, totals, status, search_method=trim(search_method))
  call require(status == 0, "false-nearest calculation failed")
  call lyapunov_stretching(series, 3, 1, 5, 2, min(500, n / 2), 10, &
    0.5_dp * sample_sd(series), stretching, refs_used, status, &
    search_method=trim(search_method), method_used=method_used)

  write(*, '(a, i0)') "observations: ", n
  write(*, '(a, a)') "neighbor search used: ", trim(method_used)
  write(*, '(a, es14.6)') "sample standard deviation: ", sample_sd(series)
  write(*, '(a, es14.6)') "AMI lag 1: ", ami(1)
  write(*, '(a, es14.6)') "C2(m=3): ", c2
  write(*, '(a, *(f8.4, 1x))') "false-neighbor fractions m=1..5: ", fractions
  if (status == 0) then
    call lyapunov_linear_fit(stretching, 2, min(7, size(stretching)), intercept, exponent, status)
    if (status == 0) then
      write(*, '(a, i0)') "Lyapunov reference points retained: ", refs_used
      write(*, '(a, es14.6)') "early stretching slope per sample: ", exponent
    end if
  else
    write(*, '(a)') "Lyapunov estimate unavailable for the selected default neighborhood."
  end if
contains
  subroutine grow(x, current_capacity)
    real(dp), allocatable, intent(inout) :: x(:)
    integer, intent(inout) :: current_capacity
    real(dp), allocatable :: tmp(:)
    allocate(tmp(2 * current_capacity))
    tmp(:current_capacity) = x
    call move_alloc(tmp, x)
    current_capacity = 2 * current_capacity
  end subroutine grow

  pure real(dp) function sample_sd(x) result(sd)
    real(dp), intent(in) :: x(:)
    real(dp) :: mean_x
    mean_x = sum(x) / real(size(x), dp)
    sd = sqrt(sum((x - mean_x)**2) / real(size(x) - 1, dp))
  end function sample_sd

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require
end program analyze_csv
