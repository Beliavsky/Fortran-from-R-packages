! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
program demo_tserieschaos
  use tserieschaos, only : dp, simulate_lorenz, average_mutual_information, &
    correlation_integral, false_nearest_curve, lyapunov_stretching, lyapunov_linear_fit
  implicit none
  real(dp), allocatable :: times(:), states(:, :), ami(:), fractions(:), stretching(:)
  integer, allocatable :: totals(:)
  real(dp) :: c2, intercept, exponent
  integer :: status, refs_used, first
  character(len=16) :: method_used

  call simulate_lorenz(0.0_dp, 30.0_dp, 0.01_dp, [1.0_dp, 1.0_dp, 1.0_dp], &
    [10.0_dp, 28.0_dp, -8.0_dp / 3.0_dp], times, states, status)
  call require(status == 0, "Lorenz simulation failed")
  first = 1001

  call average_mutual_information(states(first:, 1), 32, 20, ami, status)
  call require(status == 0, "AMI calculation failed")
  call correlation_integral(states(first:, 1), 3, 4, 20, 3.0_dp, c2, status)
  call require(status == 0, "correlation integral failed")
  call false_nearest_curve(states(first:, 1), 5, 4, 20, 10.0_dp, 8.0_dp, &
    fractions, totals, status, search_method="auto")
  call require(status == 0, "false-nearest calculation failed")
  call lyapunov_stretching(states(first:, 1), 3, 4, 20, 2, 500, 15, 8.0_dp, &
    stretching, refs_used, status, search_method="auto", method_used=method_used)
  call require(status == 0, "Lyapunov tracking failed")
  call lyapunov_linear_fit(stretching, 2, 8, intercept, exponent, status, dt=0.01_dp)
  call require(status == 0, "Lyapunov regression failed")

  write(*, '(a, i0)') "Lorenz observations: ", size(states, 1) - first + 1
  write(*, '(a, es14.6)') "AMI at lag 1: ", ami(1)
  write(*, '(a, es14.6)') "C2(m=3, eps=3): ", c2
  write(*, '(a, *(f8.4, 1x))') "False-neighbor fractions m=1..5: ", fractions
  write(*, '(a, a)') "Neighbor search used: ", trim(method_used)
  write(*, '(a, i0)') "Lyapunov reference points retained: ", refs_used
  write(*, '(a, es14.6)') "Early stretching slope per time unit: ", exponent
contains
  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require
end program demo_tserieschaos
