! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
program lorenz_analysis
  use tserieschaos, only : dp, simulate_lorenz, recurrence_distance_matrix, &
    correlation_dimension_curve
  implicit none
  real(dp), allocatable :: times(:), states(:, :), recurrence(:, :), eps(:), c2(:, :)
  integer :: status

  call simulate_lorenz(0.0_dp, 5.0_dp, 0.02_dp, [1.0_dp, 1.0_dp, 1.0_dp], &
    [10.0_dp, 28.0_dp, -8.0_dp / 3.0_dp], times, states, status)
  if (status /= 0) error stop "simulation failed"
  call recurrence_distance_matrix(states(101:, 1), 3, 2, recurrence, status)
  if (status /= 0) error stop "recurrence calculation failed"
  call correlation_dimension_curve(states(101:, 1), 4, 2, 10, 0.05_dp, 20, eps, c2, status)
  if (status /= 0) error stop "correlation-dimension calculation failed"

  write(*, '(a, i0, a, i0)') "recurrence matrix: ", size(recurrence, 1), " x ", size(recurrence, 2)
  write(*, '(a, es14.6)') "largest epsilon: ", eps(size(eps))
  write(*, '(a, es14.6)') "C2 at largest epsilon, m=4: ", c2(size(eps), 4)
end program lorenz_analysis
