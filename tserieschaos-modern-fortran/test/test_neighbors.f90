! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
program test_neighbors
  use tserieschaos, only : dp, false_nearest_fraction, false_nearest_curve, &
    find_k_nearests, follow_neighbor_points, lyapunov_stretching, &
    lyapunov_linear_fit, simulate_lorenz
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  real(dp), allocatable :: fractions(:), distances(:, :), stretching(:), times(:)
  real(dp), allocatable :: states(:, :), search_series(:)
  real(dp), allocatable :: direct_distances(:, :), box_distances(:, :)
  real(dp), allocatable :: direct_stretching(:), box_stretching(:)
  integer, allocatable :: totals(:), nearest(:, :)
  integer, allocatable :: direct_nearest(:, :), box_nearest(:, :)
  real(dp) :: fraction, direct_fraction, box_fraction, intercept, exponent
  integer :: total, direct_total, box_total, status, refs_used, box_refs, i
  integer(int64) :: direct_evaluations, box_evaluations, auto_evaluations
  real(dp) :: line(6)
  character(len=16) :: method_used

  call false_nearest_fraction(&
    [0.0_dp, 0.2_dp, 0.8_dp, 0.1_dp, 0.9_dp, 0.3_dp, 0.7_dp, 0.4_dp], &
    1, 1, 0, 10.0_dp, 0.7_dp, fraction, total, status)
  call check(status == 0, "false_nearest_fraction status")
  call check(total > 0 .and. fraction >= 0.0_dp .and. fraction <= 1.0_dp, &
    "false_nearest_fraction range")

  call false_nearest_curve(&
    [0.0_dp, 0.2_dp, 0.8_dp, 0.1_dp, 0.9_dp, 0.3_dp, 0.7_dp, 0.4_dp, &
      0.6_dp, 0.5_dp], &
    3, 1, 0, 10.0_dp, 0.8_dp, fractions, totals, status)
  call check(status == 0, "false_nearest_curve status")
  call check(size(fractions) == 3 .and. all(totals > 0), &
    "false_nearest_curve outputs")

  call find_k_nearests([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], &
    1, 1, 0, 10.0_dp, 3, 2, 1, nearest, distances, status)
  call check(status == 0, "find_k_nearests status")
  call check(all(shape(nearest) == [3, 2]), "find_k_nearests shape")
  call check(all(nearest(:, 1) /= [1, 2, 3]), "find_k_nearests excludes self")
  call check(maxval(abs(distances(:, 1) - 1.0_dp)) < 1.0e-14_dp, &
    "find_k_nearests first distance")
  call check(all(distances(:, 2) >= distances(:, 1)), "find_k_nearests sorted")

  allocate(search_series(1200))
  do i = 1, size(search_series)
    search_series(i) = real(mod(104729 * i, 1000003), dp) / 1000003.0_dp
  end do
  call find_k_nearests(search_series, 2, 1, 5, 0.08_dp, 400, 3, 2, &
    direct_nearest, direct_distances, status, search_method="direct", &
    distance_evaluations=direct_evaluations, method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "direct", &
    "direct search selection")
  call find_k_nearests(search_series, 2, 1, 5, 0.08_dp, 400, 3, 2, &
    box_nearest, box_distances, status, search_method="box", &
    distance_evaluations=box_evaluations, method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "box", "box search selection")
  call check(all(direct_nearest == box_nearest), "direct and box neighbor indices")
  call check(maxval(abs(direct_distances - box_distances)) < 1.0e-14_dp, &
    "direct and box neighbor distances")
  call check(box_evaluations < direct_evaluations / 4_int64, &
    "box search reduces exact distance evaluations")

  call find_k_nearests(search_series, 2, 1, 5, 0.08_dp, 400, 3, 2, nearest, &
    distances, status, search_method="auto", &
    distance_evaluations=auto_evaluations, method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "box", &
    "auto search chooses box for large low-dimensional data")
  call check(all(nearest == direct_nearest), "auto search neighbor equivalence")
  call check(auto_evaluations == box_evaluations, "auto search evaluation count")

  call find_k_nearests(search_series, 9, 1, 5, 0.08_dp, 20, 1, 2, nearest, &
    distances, status, search_method="auto", method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "direct", &
    "auto search falls back for high embedding dimension")

  call false_nearest_fraction(search_series, 2, 1, 5, 10.0_dp, 0.08_dp, &
    direct_fraction, direct_total, status, search_method="direct", &
    distance_evaluations=direct_evaluations)
  call check(status == 0, "direct false-neighbor search")
  call false_nearest_fraction(search_series, 2, 1, 5, 10.0_dp, 0.08_dp, &
    box_fraction, box_total, status, search_method="box", &
    distance_evaluations=box_evaluations)
  call check(status == 0, "box false-neighbor search")
  call check(direct_total == box_total, "direct and box false-neighbor totals")
  call check(abs(direct_fraction - box_fraction) < 1.0e-15_dp, &
    "direct and box false-neighbor fractions")
  call check(box_evaluations < direct_evaluations / 4_int64, &
    "box false-neighbor search reduces evaluations")

  deallocate(nearest)
  allocate(nearest(1, 1))
  nearest(1, 1) = 2
  call follow_neighbor_points([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
    1, 1, [1], nearest, 3, stretching, status)
  call check(status == 0, "follow_neighbor_points status")
  call check(maxval(abs(stretching)) < 1.0e-14_dp, "follow_neighbor_points exact")

  do i = 1, 6
    line(i) = 2.0_dp + 0.5_dp * real(i - 1, dp)
  end do
  call lyapunov_linear_fit(line, 1, 6, intercept, exponent, status)
  call check(status == 0, "lyapunov_linear_fit status")
  call check(abs(intercept - 2.0_dp) < 1.0e-14_dp, "lyapunov intercept")
  call check(abs(exponent - 0.5_dp) < 1.0e-14_dp, "lyapunov exponent")

  call simulate_lorenz(0.0_dp, 20.0_dp, 0.01_dp, [1.0_dp, 1.0_dp, 1.0_dp], &
    [10.0_dp, 28.0_dp, -8.0_dp / 3.0_dp], times, states, status)
  call check(status == 0, "Lorenz simulation for Lyapunov")
  call lyapunov_stretching(states(:, 1), 3, 4, 20, 2, 500, 12, 8.0_dp, &
    direct_stretching, refs_used, status, search_method="direct", &
    distance_evaluations=direct_evaluations)
  call check(status == 0, "direct lyapunov_stretching status")
  call lyapunov_stretching(states(:, 1), 3, 4, 20, 2, 500, 12, 8.0_dp, &
    box_stretching, box_refs, status, search_method="box", &
    distance_evaluations=box_evaluations, method_used=method_used)
  call check(status == 0 .and. trim(method_used) == "box", &
    "box lyapunov_stretching status")
  call check(refs_used == box_refs, "direct and box Lyapunov reference counts")
  call check(maxval(abs(direct_stretching - box_stretching)) < 1.0e-14_dp, &
    "direct and box Lyapunov paths")
  call check(box_evaluations < direct_evaluations, &
    "box Lyapunov search reduces evaluations")
  call check(refs_used > 20 .and. all(ieee_is_finite(box_stretching)), &
    "lyapunov_stretching outputs")

  print '(a)', "Direct, box-index, auto-search, and Lyapunov tests passed."
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
