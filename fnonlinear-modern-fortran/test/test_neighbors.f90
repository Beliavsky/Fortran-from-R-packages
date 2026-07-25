! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program test_neighbors
  use fnonlinear, only : dp, false_nearest_neighbors, find_k_nearests, &
    lyapunov_stretching, lyapunov_linear_fit, lorenz_sim
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  real(dp), allocatable :: fractions(:), distances(:, :), search_series(:)
  real(dp), allocatable :: direct_distances(:, :), box_distances(:, :)
  real(dp), allocatable :: direct_stretching(:), box_stretching(:)
  real(dp), allocatable :: times(:), trajectory(:, :)
  integer, allocatable :: totals(:), nearest(:, :), direct_nearest(:, :), box_nearest(:, :)
  real(dp) :: intercept, exponent
  integer :: status, refs_used, box_refs, i
  integer(int64) :: direct_evaluations, box_evaluations
  character(len=16) :: method_used

  call false_nearest_neighbors([0.0_dp, 0.2_dp, 0.8_dp, 0.1_dp, 0.9_dp, &
    0.3_dp, 0.7_dp, 0.4_dp, 0.6_dp, 0.5_dp], 3, 1, 0, 10.0_dp, 0.8_dp, &
    fractions, totals, status, search_method="direct")
  call check(status == 0 .and. size(fractions) == 3 .and. all(totals > 0), &
    "false-nearest curve")
  call check(all(fractions >= 0.0_dp .and. fractions <= 1.0_dp), &
    "false-nearest range")

  allocate(search_series(1200))
  do i = 1, size(search_series)
    search_series(i) = real(mod(104729 * i, 1000003), dp) / 1000003.0_dp
  end do
  call find_k_nearests(search_series, 2, 1, 5, 0.08_dp, 400, 3, 2, &
    direct_nearest, direct_distances, status, search_method="direct", &
    distance_evaluations=direct_evaluations, method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "direct", "direct search")
  call find_k_nearests(search_series, 2, 1, 5, 0.08_dp, 400, 3, 2, &
    box_nearest, box_distances, status, search_method="box", &
    distance_evaluations=box_evaluations, method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "box", "box search")
  call check(all(direct_nearest == box_nearest), "neighbor-index equivalence")
  call check(maxval(abs(direct_distances - box_distances)) < 1.0e-14_dp, &
    "neighbor-distance equivalence")
  call check(box_evaluations < direct_evaluations, "box acceleration")

  call find_k_nearests(search_series, 2, 1, 5, 0.08_dp, 50, 2, 2, nearest, &
    distances, status, search_method="auto", method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "box", "automatic box search")

  call lyapunov_linear_fit([2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp, 4.0_dp, 4.5_dp], &
    1, 6, intercept, exponent, status)
  call check(status == 0 .and. abs(intercept - 2.0_dp) < 1.0e-14_dp, &
    "Lyapunov intercept")
  call check(abs(exponent - 0.5_dp) < 1.0e-14_dp, "Lyapunov slope")

  allocate(times(2001))
  do i = 1, size(times)
    times(i) = 0.01_dp * real(i - 1, dp)
  end do
  call lorenz_sim(times, [10.0_dp, 28.0_dp, 8.0_dp / 3.0_dp], &
    [1.0_dp, 1.0_dp, 1.0_dp], trajectory, status)
  call check(status == 0, "Lorenz series")
  call lyapunov_stretching(trajectory(:, 2), 3, 4, 20, 2, 500, 12, 8.0_dp, &
    direct_stretching, refs_used, status, search_method="direct", &
    distance_evaluations=direct_evaluations)
  call check(status == 0, "direct Lyapunov path")
  call lyapunov_stretching(trajectory(:, 2), 3, 4, 20, 2, 500, 12, 8.0_dp, &
    box_stretching, box_refs, status, search_method="box", &
    distance_evaluations=box_evaluations)
  call check(status == 0, "box Lyapunov path")
  call check(refs_used == box_refs .and. refs_used > 20, "Lyapunov reference count")
  call check(maxval(abs(direct_stretching - box_stretching)) < 1.0e-14_dp, &
    "Lyapunov method equivalence")
  call check(all(ieee_is_finite(box_stretching)), "finite Lyapunov path")

  print '(a)', "False-neighbor, box-index, and Lyapunov tests passed."
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // message
      error stop 1
    end if
  end subroutine check
end program test_neighbors
