! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1988-1990 Blake LeBaron
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program analyze_csv
  use fnonlinear, only : dp, rng_state, rng_seed, mutual_information_curve, &
    correlation_integral, false_nearest_neighbors, lyapunov_stretching, &
    lyapunov_linear_fit, bds_test_result, neural_test_result, runs_test_result, &
    bds_test, white_neural_test, terasvirta_neural_test, runs_test
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  character(len=512) :: filename, line
  character(len=16) :: search_method, method_used
  real(dp), allocatable :: series(:), ami(:), fractions(:), stretching(:)
  integer, allocatable :: totals(:)
  type(bds_test_result) :: bds
  type(neural_test_result) :: white, tera
  type(runs_test_result) :: runs
  type(rng_state) :: rng
  real(dp) :: c2, intercept, exponent, value, sd, centered_mean
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
  if (n < 80) error stop "at least 80 observations are required"
  series = series(:n)
  sd = sample_sd(series)
  centered_mean = sum(series) / real(n, dp)

  call mutual_information_curve(series, 16, 20, ami, status)
  call require(status == 0, "mutual-information calculation failed")
  call correlation_integral(series, 3, 1, 5, 0.1_dp * sd, c2, status)
  call require(status == 0, "correlation-integral calculation failed")
  call false_nearest_neighbors(series, 5, 1, 5, 10.0_dp, 0.2_dp * sd, &
    fractions, totals, status, search_method=trim(search_method))
  call require(status == 0, "false-nearest calculation failed")
  call bds_test(series, 3, [0.5_dp * sd, sd, 1.5_dp * sd, 2.0_dp * sd], bds, status)
  call require(status == 0, "BDS test failed")
  call rng_seed(rng, 20260725_int64)
  call white_neural_test(series, 2, 2, 10, 4.0_dp, white, status, rng=rng)
  call require(status == 0, "White neural-network test failed")
  call terasvirta_neural_test(series, 2, tera, status)
  call require(status == 0, "Terasvirta neural-network test failed")
  call runs_test(series - centered_mean, runs, status)
  call require(status == 0, "runs test failed")
  call lyapunov_stretching(series, 3, 1, 5, 2, min(500, n / 2), 10, &
    0.5_dp * sd, stretching, refs_used, status, search_method=trim(search_method), &
    method_used=method_used)

  write(*, '(a, i0)') "observations: ", n
  write(*, '(a, es14.6)') "sample standard deviation: ", sd
  write(*, '(a, es14.6)') "mutual information lag 1: ", ami(1)
  write(*, '(a, es14.6)') "correlation integral C2(m=3): ", c2
  write(*, '(a, *(f8.4, 1x))') "false-neighbor fractions m=1..5: ", fractions
  write(*, '(a, *(f10.4, 1x))') "BDS statistics m=2..3 at epsilon=sd: ", bds%statistic(:, 2)
  write(*, '(a, es14.6)') "White chi-square p-value: ", white%chi_square_p
  write(*, '(a, es14.6)') "Terasvirta chi-square p-value: ", tera%chi_square_p
  write(*, '(a, es14.6)') "runs-test p-value: ", runs%p_value
  if (status == 0) then
    call lyapunov_linear_fit(stretching, 2, min(7, size(stretching)), intercept, exponent, status)
    if (status == 0) then
      write(*, '(a, a)') "neighbor search used: ", trim(method_used)
      write(*, '(a, i0)') "Lyapunov reference points retained: ", refs_used
      write(*, '(a, es14.6)') "early stretching slope per sample: ", exponent
    end if
  else
    write(*, '(a)') "Lyapunov estimate unavailable for the default neighborhood."
  end if
contains
  subroutine grow(x, current_capacity)
    real(dp), allocatable, intent(inout) :: x(:)
    integer, intent(inout) :: current_capacity
    real(dp), allocatable :: temporary(:)
    allocate(temporary(2 * current_capacity))
    temporary(:current_capacity) = x
    call move_alloc(temporary, x)
    current_capacity = 2 * current_capacity
  end subroutine grow

  pure real(dp) function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: mean_value
    mean_value = sum(x) / real(size(x), dp)
    value = sqrt(sum((x - mean_value)**2) / real(size(x) - 1, dp))
  end function sample_sd

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require
end program analyze_csv
