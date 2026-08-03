! SPDX-License-Identifier: GPL-3.0-only
module mixedind_rng
  use iso_fortran_env, only : int64
  use mixedind_kinds, only : dp
  implicit none
  private

  type, public :: rng_state
    integer(int64) :: state = 88172645463393265_int64
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

  public :: rng_seed, rng_uniform, rng_normal, rng_gamma, rng_poisson

contains

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(int64), intent(in) :: seed

    if (seed == 0_int64) then
      rng%state = 88172645463393265_int64
    else
      rng%state = seed
    end if
    rng%has_spare = .false.
    rng%spare = 0.0_dp
  end subroutine rng_seed

  function next_u64(rng) result(x)
    type(rng_state), intent(inout) :: rng
    integer(int64) :: x

    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
  end function next_u64

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(int64) :: x, mask

    mask = int(z'001FFFFFFFFFFFFF', int64)
    x = iand(next_u64(rng), mask)
    u = (real(x, dp) + 0.5_dp) / 9007199254740992.0_dp
    u = min(max(u, epsilon(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
  end function rng_uniform

  function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z
    real(dp) :: u1, u2, r, theta
    real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp

    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if

    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    r = sqrt(-2.0_dp * log(u1))
    theta = two_pi * u2
    z = r * cos(theta)
    rng%spare = r * sin(theta)
    rng%has_spare = .true.
  end function rng_normal

  recursive function rng_gamma(rng, shape) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp) :: x
    real(dp) :: d, c, z, v, u

    if (shape <= 0.0_dp) then
      x = 0.0_dp
      return
    end if

    if (shape < 1.0_dp) then
      x = rng_gamma(rng, shape + 1.0_dp) * rng_uniform(rng)**(1.0_dp / shape)
      return
    end if

    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      z = rng_normal(rng)
      v = 1.0_dp + c * z
      if (v <= 0.0_dp) cycle
      v = v * v * v
      u = rng_uniform(rng)
      if (u < 1.0_dp - 0.0331_dp * z**4) exit
      if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
    end do
    x = d * v
  end function rng_gamma

  function rng_poisson(rng, lambda) result(k)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: lambda
    integer :: k
    real(dp) :: p, l, z, x

    if (lambda <= 0.0_dp) then
      k = 0
    else if (lambda < 30.0_dp) then
      l = exp(-lambda)
      p = 1.0_dp
      k = 0
      do
        k = k + 1
        p = p * rng_uniform(rng)
        if (p <= l) exit
      end do
      k = k - 1
    else
      do
        z = rng_normal(rng)
        x = lambda + sqrt(lambda) * z
        if (x >= 0.0_dp) exit
      end do
      k = nint(x)
    end if
  end function rng_poisson

end module mixedind_rng
