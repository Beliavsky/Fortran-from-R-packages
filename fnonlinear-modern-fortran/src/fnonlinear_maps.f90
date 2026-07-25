! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear_maps
  use chaos_kinds, only : dp
  use fnonlinear_rng, only : rng_state, rng_uniform
  implicit none
  private

  abstract interface
    subroutine ode_rhs(t, state, parms, deriv)
      import dp
      real(dp), intent(in) :: t, state(:), parms(:)
      real(dp), intent(out) :: deriv(:)
    end subroutine ode_rhs
  end interface

  public :: ode_rhs, rk4_integrate_times
  public :: tent_sim, henon_sim, ikeda_sim, logistic_sim
  public :: lorenz_rhs, rossler_rhs, lorenz_sim, rossler_sim
  public :: lorentz_sim, roessler_sim
contains
  subroutine tent_sim(n, n_skip, a, x, status, start, rng)
    integer, intent(in) :: n, n_skip
    real(dp), intent(in) :: a
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: start
    type(rng_state), intent(inout), optional :: rng
    real(dp), allocatable :: work(:)
    real(dp) :: x0, aa
    type(rng_state) :: local_rng
    integer :: i, ntotal

    if (n < 1 .or. n_skip < 0 .or. a <= 0.0_dp) then
      allocate(x(0))
      status = 1
      return
    end if
    if (present(start)) then
      x0 = start
    else if (present(rng)) then
      x0 = rng_uniform(rng)
    else
      x0 = rng_uniform(local_rng)
    end if
    aa = a
    if (abs(aa - 2.0_dp) <= epsilon(aa)) aa = aa - epsilon(aa)
    ntotal = n + n_skip
    allocate(work(ntotal), x(n))
    work(1) = x0
    do i = 2, ntotal
      work(i) = 0.5_dp * aa * (1.0_dp - 2.0_dp * abs(work(i - 1) - 0.5_dp))
    end do
    x = work(n_skip + 1:ntotal)
    status = 0
  end subroutine tent_sim

  subroutine logistic_sim(n, n_skip, r, x, status, start, rng)
    integer, intent(in) :: n, n_skip
    real(dp), intent(in) :: r
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: start
    type(rng_state), intent(inout), optional :: rng
    real(dp), allocatable :: work(:)
    real(dp) :: x0
    type(rng_state) :: local_rng
    integer :: i, ntotal

    if (n < 1 .or. n_skip < 0 .or. r <= 0.0_dp) then
      allocate(x(0))
      status = 1
      return
    end if
    if (present(start)) then
      x0 = start
    else if (present(rng)) then
      x0 = rng_uniform(rng)
    else
      x0 = rng_uniform(local_rng)
    end if
    ntotal = n + n_skip
    allocate(work(ntotal), x(n))
    work(1) = x0
    do i = 2, ntotal
      work(i) = r * work(i - 1) * (1.0_dp - work(i - 1))
    end do
    x = work(n_skip + 1:ntotal)
    status = 0
  end subroutine logistic_sim

  subroutine henon_sim(n, n_skip, a, b, xy, status, start, rng)
    integer, intent(in) :: n, n_skip
    real(dp), intent(in) :: a, b
    real(dp), allocatable, intent(out) :: xy(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: start(:)
    type(rng_state), intent(inout), optional :: rng
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: x0, y0
    type(rng_state) :: local_rng
    integer :: i, ntotal

    if (n < 1 .or. n_skip < 0) then
      allocate(xy(0, 0))
      status = 1
      return
    end if
    if (present(start)) then
      if (size(start) < 2) then
        allocate(xy(0, 0))
        status = 2
        return
      end if
      x0 = start(1)
      y0 = start(2)
    else if (present(rng)) then
      x0 = rng_uniform(rng)
      y0 = rng_uniform(rng)
    else
      x0 = rng_uniform(local_rng)
      y0 = rng_uniform(local_rng)
    end if
    ntotal = n + n_skip
    allocate(x(ntotal), y(ntotal), xy(n, 2))
    x(1) = x0
    y(1) = y0
    do i = 2, ntotal
      x(i) = 1.0_dp - a * x(i - 1)**2 + b * y(i - 1)
      y(i) = x(i - 1)
    end do
    xy(:, 1) = x(n_skip + 1:ntotal)
    xy(:, 2) = y(n_skip + 1:ntotal)
    status = 0
  end subroutine henon_sim

  subroutine ikeda_sim(n, n_skip, a, b, c, z_out, status, start, rng)
    integer, intent(in) :: n, n_skip
    real(dp), intent(in) :: a, b, c
    real(dp), allocatable, intent(out) :: z_out(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: start(:)
    type(rng_state), intent(inout), optional :: rng
    complex(dp), allocatable :: z(:)
    complex(dp) :: z0, ia, ib
    type(rng_state) :: local_rng
    integer :: i, ntotal

    if (n < 1 .or. n_skip < 0) then
      allocate(z_out(0, 0))
      status = 1
      return
    end if
    if (present(start)) then
      if (size(start) < 2) then
        allocate(z_out(0, 0))
        status = 2
        return
      end if
      z0 = cmplx(start(1), start(2), kind=dp)
    else if (present(rng)) then
      z0 = cmplx(rng_uniform(rng), rng_uniform(rng), kind=dp)
    else
      z0 = cmplx(rng_uniform(local_rng), rng_uniform(local_rng), kind=dp)
    end if
    ia = cmplx(0.0_dp, a, kind=dp)
    ib = cmplx(0.0_dp, b, kind=dp)
    ntotal = n + n_skip
    allocate(z(ntotal), z_out(n, 2))
    z(1) = z0
    do i = 2, ntotal
      z(i) = 1.0_dp + c * z(i - 1) * exp(ia - ib / (1.0_dp + abs(z(i - 1))**2))
    end do
    z_out(:, 1) = real(z(n_skip + 1:ntotal), dp)
    z_out(:, 2) = aimag(z(n_skip + 1:ntotal))
    status = 0
  end subroutine ikeda_sim

  subroutine rk4_integrate_times(rhs, start_state, times, parms, states, status)
    procedure(ode_rhs) :: rhs
    real(dp), intent(in) :: start_state(:), times(:), parms(:)
    real(dp), allocatable, intent(out) :: states(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), temp(:)
    real(dp) :: h, t
    integer :: i, nstate

    if (size(start_state) < 1 .or. size(times) < 1) then
      allocate(states(0, 0))
      status = 1
      return
    end if
    if (size(times) > 1) then
      if (any(times(2:) <= times(:size(times) - 1))) then
        allocate(states(0, 0))
        status = 2
        return
      end if
    end if
    nstate = size(start_state)
    allocate(states(size(times), nstate))
    allocate(k1(nstate), k2(nstate), k3(nstate), k4(nstate), temp(nstate))
    states(1, :) = start_state
    do i = 1, size(times) - 1
      t = times(i)
      h = times(i + 1) - times(i)
      call rhs(t, states(i, :), parms, k1)
      temp = states(i, :) + 0.5_dp * h * k1
      call rhs(t + 0.5_dp * h, temp, parms, k2)
      temp = states(i, :) + 0.5_dp * h * k2
      call rhs(t + 0.5_dp * h, temp, parms, k3)
      temp = states(i, :) + h * k3
      call rhs(t + h, temp, parms, k4)
      states(i + 1, :) = states(i, :) + h * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4) / 6.0_dp
    end do
    status = 0
  end subroutine rk4_integrate_times

  subroutine lorenz_rhs(t, state, parms, deriv)
    real(dp), intent(in) :: t, state(:), parms(:)
    real(dp), intent(out) :: deriv(:)
    if (size(state) /= 3 .or. size(parms) < 3 .or. size(deriv) /= 3) then
      error stop "lorenz_rhs: invalid dimensions"
    end if
    deriv(1) = parms(1) * (state(2) - state(1)) + 0.0_dp * t
    deriv(2) = -state(1) * state(3) + parms(2) * state(1) - state(2)
    deriv(3) = state(1) * state(2) - parms(3) * state(3)
  end subroutine lorenz_rhs

  subroutine rossler_rhs(t, state, parms, deriv)
    real(dp), intent(in) :: t, state(:), parms(:)
    real(dp), intent(out) :: deriv(:)
    if (size(state) /= 3 .or. size(parms) < 3 .or. size(deriv) /= 3) then
      error stop "rossler_rhs: invalid dimensions"
    end if
    deriv(1) = -(state(2) + state(3)) + 0.0_dp * t
    deriv(2) = state(1) + parms(1) * state(2)
    deriv(3) = parms(2) + state(1) * state(3) - parms(3) * state(3)
  end subroutine rossler_rhs

  subroutine lorenz_sim(times, parms, start, trajectory, status)
    real(dp), intent(in) :: times(:), parms(:), start(:)
    real(dp), allocatable, intent(out) :: trajectory(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: states(:, :)
    call rk4_integrate_times(lorenz_rhs, start, times, parms, states, status)
    if (status /= 0) then
      allocate(trajectory(0, 0))
      return
    end if
    allocate(trajectory(size(times), 4))
    trajectory(:, 1) = times
    trajectory(:, 2:4) = states
  end subroutine lorenz_sim

  subroutine rossler_sim(times, parms, start, trajectory, status)
    real(dp), intent(in) :: times(:), parms(:), start(:)
    real(dp), allocatable, intent(out) :: trajectory(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: states(:, :)
    call rk4_integrate_times(rossler_rhs, start, times, parms, states, status)
    if (status /= 0) then
      allocate(trajectory(0, 0))
      return
    end if
    allocate(trajectory(size(times), 4))
    trajectory(:, 1) = times
    trajectory(:, 2:4) = states
  end subroutine rossler_sim

  subroutine lorentz_sim(times, parms, start, trajectory, status)
    real(dp), intent(in) :: times(:), parms(:), start(:)
    real(dp), allocatable, intent(out) :: trajectory(:, :)
    integer, intent(out) :: status
    call lorenz_sim(times, parms, start, trajectory, status)
  end subroutine lorentz_sim

  subroutine roessler_sim(times, parms, start, trajectory, status)
    real(dp), intent(in) :: times(:), parms(:), start(:)
    real(dp), allocatable, intent(out) :: trajectory(:, :)
    integer, intent(out) :: status
    call rossler_sim(times, parms, start, trajectory, status)
  end subroutine roessler_sim
end module fnonlinear_maps
