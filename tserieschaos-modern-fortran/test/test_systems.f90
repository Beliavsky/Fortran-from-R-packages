! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
module test_observation_mod
  use chaos_kinds, only : dp
  implicit none
contains
  function phase_observation(state) result(value)
    real(dp), intent(in) :: state(:)
    real(dp) :: value
    value = 2.0_dp * state(3)
  end function phase_observation
end module test_observation_mod

program test_systems
  use tserieschaos, only : dp, lorenz_rhs, rossler_rhs, duffing_rhs, simulate_observed, &
    simulate_lorenz, simulate_rossler, simulate_duffing
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use test_observation_mod, only : phase_observation
  implicit none
  real(dp) :: deriv(3)
  real(dp), allocatable :: times(:), states(:, :), observed(:)
  integer :: status

  call lorenz_rhs(0.0_dp, [1.0_dp, 2.0_dp, 3.0_dp], [10.0_dp, 28.0_dp, -8.0_dp / 3.0_dp], deriv)
  call check(maxval(abs(deriv - [10.0_dp, 23.0_dp, -6.0_dp])) < 1.0e-14_dp, "Lorenz RHS")

  call rossler_rhs(0.0_dp, [1.0_dp, 2.0_dp, 3.0_dp], [0.2_dp, 0.2_dp, 5.7_dp], deriv)
  call check(maxval(abs(deriv - [-5.0_dp, 1.4_dp, -13.9_dp])) < 1.0e-14_dp, "Rossler RHS")

  call duffing_rhs(0.0_dp, [1.0_dp, 2.0_dp, 0.0_dp], [0.3_dp, 0.5_dp, 1.2_dp], deriv)
  call check(maxval(abs(deriv - [2.0_dp, -0.1_dp, 1.2_dp])) < 1.0e-14_dp, "Duffing RHS")

  call simulate_lorenz(0.0_dp, 1.0_dp, 0.01_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    [10.0_dp, 28.0_dp, -8.0_dp / 3.0_dp], times, states, status)
  call check(status == 0 .and. size(times) == 101, "Lorenz integration shape")
  call check(maxval(abs(states)) < 1.0e-14_dp, "Lorenz equilibrium")

  call simulate_duffing(0.0_dp, 2.0_dp, 0.01_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    [0.3_dp, 0.0_dp, 1.5_dp], times, states, status)
  call check(status == 0, "Duffing integration status")
  call check(maxval(abs(states(:, 1:2))) < 1.0e-13_dp, "Duffing zero orbit")
  call check(abs(states(size(states, 1), 3) - 3.0_dp) < 1.0e-12_dp, "Duffing phase")

  call simulate_rossler(0.0_dp, 5.0_dp, 0.01_dp, [1.0_dp, 0.0_dp, 0.0_dp], &
    [0.2_dp, 0.2_dp, 5.7_dp], times, states, status)
  call check(status == 0 .and. all(ieee_is_finite(states)), "Rossler finite integration")

  call simulate_observed(duffing_rhs, 0.0_dp, 1.0_dp, 0.01_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    [0.3_dp, 0.0_dp, 2.0_dp], times, observed, status, phase_observation)
  call check(status == 0, "observed simulation status")
  call check(abs(observed(size(observed)) - 4.0_dp) < 1.0e-12_dp, "custom observation callback")

  print '(a)', "Continuous-system and RK4 tests passed."
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // message
      error stop 1
    end if
  end subroutine check
end program test_systems
