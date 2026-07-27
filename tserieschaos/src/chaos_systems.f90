! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
module chaos_systems
  use chaos_kinds, only : dp
  implicit none
  private
  public :: rhs_proc, observation_proc
  public :: lorenz_rhs, rossler_rhs, duffing_rhs, integrate_rk4, simulate_observed
  public :: simulate_lorenz, simulate_rossler, simulate_duffing

  abstract interface
    subroutine rhs_proc(t, state, parms, deriv)
      import dp
      real(dp), intent(in) :: t
      real(dp), intent(in) :: state(:), parms(:)
      real(dp), intent(out) :: deriv(:)
    end subroutine rhs_proc

    function observation_proc(state) result(value)
      import dp
      real(dp), intent(in) :: state(:)
      real(dp) :: value
    end function observation_proc
  end interface
contains
  subroutine lorenz_rhs(t, state, parms, deriv)
    real(dp), intent(in) :: t, state(:), parms(:)
    real(dp), intent(out) :: deriv(:)
    if (size(state) /= 3 .or. size(parms) < 3 .or. size(deriv) /= 3) error stop "lorenz_rhs: invalid size"
    deriv(1) = parms(1) * (state(2) - state(1)) + 0.0_dp * t
    deriv(2) = state(1) * (parms(2) - state(3)) - state(2)
    deriv(3) = state(1) * state(2) + parms(3) * state(3)
  end subroutine lorenz_rhs

  subroutine rossler_rhs(t, state, parms, deriv)
    real(dp), intent(in) :: t, state(:), parms(:)
    real(dp), intent(out) :: deriv(:)
    if (size(state) /= 3 .or. size(parms) < 3 .or. size(deriv) /= 3) error stop "rossler_rhs: invalid size"
    deriv(1) = -(state(2) + state(3)) + 0.0_dp * t
    deriv(2) = state(1) + parms(1) * state(2)
    deriv(3) = parms(2) + state(3) * (state(1) - parms(3))
  end subroutine rossler_rhs

  subroutine duffing_rhs(t, state, parms, deriv)
    real(dp), intent(in) :: t, state(:), parms(:)
    real(dp), intent(out) :: deriv(:)
    if (size(state) /= 3 .or. size(parms) < 3 .or. size(deriv) /= 3) error stop "duffing_rhs: invalid size"
    deriv(1) = state(2) + 0.0_dp * t
    deriv(2) = state(1) * (1.0_dp - state(1)**2) - parms(1) * state(2) + parms(2) * cos(state(3))
    deriv(3) = parms(3)
  end subroutine duffing_rhs

  subroutine integrate_rk4(rhs, start_time, end_time, dt, start_state, parms, times, states, status)
    procedure(rhs_proc) :: rhs
    real(dp), intent(in) :: start_time, end_time, dt
    real(dp), intent(in) :: start_state(:), parms(:)
    real(dp), allocatable, intent(out) :: times(:), states(:, :)
    integer, intent(out) :: status
    integer :: nstep, i, nstate
    real(dp) :: t, h
    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), tmp(:)

    if (dt <= 0.0_dp .or. end_time < start_time .or. size(start_state) == 0) then
      allocate(times(0), states(0, 0))
      status = 1
      return
    end if
    nstep = int(floor((end_time - start_time) / dt + 1.0e-12_dp))
    nstate = size(start_state)
    allocate(times(nstep + 1), states(nstep + 1, nstate))
    allocate(k1(nstate), k2(nstate), k3(nstate), k4(nstate), tmp(nstate))
    times(1) = start_time
    states(1, :) = start_state
    do i = 1, nstep
      t = times(i)
      h = min(dt, end_time - t)
      call rhs(t, states(i, :), parms, k1)
      tmp = states(i, :) + 0.5_dp * h * k1
      call rhs(t + 0.5_dp * h, tmp, parms, k2)
      tmp = states(i, :) + 0.5_dp * h * k2
      call rhs(t + 0.5_dp * h, tmp, parms, k3)
      tmp = states(i, :) + h * k3
      call rhs(t + h, tmp, parms, k4)
      states(i + 1, :) = states(i, :) + h * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
      times(i + 1) = t + h
    end do
    status = 0
  end subroutine integrate_rk4


  subroutine simulate_observed(rhs, start_time, end_time, dt, start_state, parms, times, observed, status, observe)
    procedure(rhs_proc) :: rhs
    real(dp), intent(in) :: start_time, end_time, dt, start_state(:), parms(:)
    real(dp), allocatable, intent(out) :: times(:), observed(:)
    integer, intent(out) :: status
    procedure(observation_proc), optional :: observe
    real(dp), allocatable :: states(:, :)
    integer :: i

    call integrate_rk4(rhs, start_time, end_time, dt, start_state, parms, times, states, status)
    if (status /= 0) then
      allocate(observed(0))
      return
    end if
    allocate(observed(size(times)))
    if (present(observe)) then
      do i = 1, size(times)
        observed(i) = observe(states(i, :))
      end do
    else
      observed = states(:, 1)
    end if
  end subroutine simulate_observed

  subroutine simulate_lorenz(start_time, end_time, dt, start_state, parms, times, states, status)
    real(dp), intent(in) :: start_time, end_time, dt, start_state(:), parms(:)
    real(dp), allocatable, intent(out) :: times(:), states(:, :)
    integer, intent(out) :: status
    call integrate_rk4(lorenz_rhs, start_time, end_time, dt, start_state, parms, times, states, status)
  end subroutine simulate_lorenz

  subroutine simulate_rossler(start_time, end_time, dt, start_state, parms, times, states, status)
    real(dp), intent(in) :: start_time, end_time, dt, start_state(:), parms(:)
    real(dp), allocatable, intent(out) :: times(:), states(:, :)
    integer, intent(out) :: status
    call integrate_rk4(rossler_rhs, start_time, end_time, dt, start_state, parms, times, states, status)
  end subroutine simulate_rossler

  subroutine simulate_duffing(start_time, end_time, dt, start_state, parms, times, states, status)
    real(dp), intent(in) :: start_time, end_time, dt, start_state(:), parms(:)
    real(dp), allocatable, intent(out) :: times(:), states(:, :)
    integer, intent(out) :: status
    call integrate_rk4(duffing_rhs, start_time, end_time, dt, start_state, parms, times, states, status)
  end subroutine simulate_duffing
end module chaos_systems
