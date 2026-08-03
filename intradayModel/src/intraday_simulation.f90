! SPDX-License-Identifier: Apache-2.0
module intraday_simulation
  use, intrinsic :: iso_fortran_env, only : int64
  use intraday_kinds, only : dp
  use intraday_types, only : volume_parameters, intraday_ok, intraday_invalid_input
  use intraday_utils, only : parameters_valid
  implicit none
  private

  public :: simulate_intraday_volume

contains

  subroutine simulate_intraday_volume(par, n_day, volume, states, seed, status)
    type(volume_parameters), intent(in) :: par
    integer, intent(in) :: n_day
    real(dp), allocatable, intent(out) :: volume(:, :)
    real(dp), allocatable, intent(out), optional :: states(:, :)
    integer, intent(in), optional :: seed
    integer, intent(out), optional :: status

    integer :: n_bin, day, bin, t, n
    integer(int64) :: rng_state
    real(dp) :: x(2), z1, z2, observation

    if (present(status)) status = intraday_ok
    if (.not. allocated(par%phi)) then
      allocate(volume(0, 0))
      if (present(states)) allocate(states(0, 0))
      if (present(status)) status = intraday_invalid_input
      return
    end if
    n_bin = size(par%phi)
    if (n_day < 1 .or. .not. parameters_valid(par, n_bin)) then
      allocate(volume(0, 0))
      if (present(states)) allocate(states(0, 0))
      if (present(status)) status = intraday_invalid_input
      return
    end if

    rng_state = 13579_int64
    if (present(seed)) rng_state = max(1_int64, modulo(int(seed, int64), 2147483646_int64))
    allocate(volume(n_bin, n_day))
    n = n_bin * n_day
    if (present(states)) allocate(states(2, n))
    x = par%x0
    t = 0

    do day = 1, n_day
      do bin = 1, n_bin
        t = t + 1
        if (t > 1) then
          call normal_pair(rng_state, z1, z2)
          if (bin == 1) then
            x(1) = par%a_eta * x(1) + sqrt(max(0.0_dp, par%var_eta)) * z1
          end if
          x(2) = par%a_mu * x(2) + sqrt(max(0.0_dp, par%var_mu)) * z2
        end if
        call normal_pair(rng_state, z1, z2)
        observation = sum(x) + par%phi(bin) + sqrt(max(0.0_dp, par%r)) * z1
        volume(bin, day) = exp(observation)
        if (present(states)) states(:, t) = x
      end do
    end do
  end subroutine simulate_intraday_volume

  subroutine normal_pair(state, z1, z2)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: z1, z2
    real(dp) :: u1, u2, radius, angle

    call uniform01(state, u1)
    call uniform01(state, u2)
    u1 = max(u1, 1.0e-14_dp)
    radius = sqrt(-2.0_dp * log(u1))
    angle = 2.0_dp * acos(-1.0_dp) * u2
    z1 = radius * cos(angle)
    z2 = radius * sin(angle)
  end subroutine normal_pair

  subroutine uniform01(state, value)
    integer(int64), intent(inout) :: state
    real(dp), intent(out) :: value
    integer(int64), parameter :: modulus = 2147483647_int64
    integer(int64), parameter :: multiplier = 16807_int64

    state = modulo(multiplier * state, modulus)
    if (state <= 0_int64) state = 1_int64
    value = real(state, dp) / real(modulus, dp)
  end subroutine uniform01

end module intraday_simulation
