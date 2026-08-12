! cmaes-fortran - GPL-2.0-only
module cmaes_rng
  use iso_fortran_env, only : int64
  use cmaes_kinds, only : dp
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal

  type :: rng_state
    integer(int64) :: state = 1_int64
    logical :: have_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(int64), intent(in) :: seed
    integer(int64), parameter :: m = 2147483647_int64
    rng%state = modulo(abs(seed), m - 1_int64) + 1_int64
    rng%have_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(int64), parameter :: a = 16807_int64
    integer(int64), parameter :: m = 2147483647_int64
    rng%state = modulo(a * rng%state, m)
    if (rng%state <= 0_int64) rng%state = 1_int64
    u = real(rng%state, dp) / real(m, dp)
  end function rng_uniform

  function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z
    real(dp) :: u1, u2, r, theta
    real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)

    if (rng%have_spare) then
      z = rng%spare
      rng%have_spare = .false.
      return
    end if

    u1 = max(rng_uniform(rng), tiny(1.0_dp))
    u2 = rng_uniform(rng)
    r = sqrt(-2.0_dp * log(u1))
    theta = twopi * u2
    z = r * cos(theta)
    rng%spare = r * sin(theta)
    rng%have_spare = .true.
  end function rng_normal
end module cmaes_rng
