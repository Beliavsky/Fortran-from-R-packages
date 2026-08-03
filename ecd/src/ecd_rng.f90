! SPDX-License-Identifier: Artistic-2.0
module ecd_rng
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use ecd_kinds, only : dp, i8, pi
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_exponential, rng_gamma

  type :: rng_state
    integer(i8) :: state = 88172645463393265_i8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

contains

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed
    if (seed == 0_i8) then
      rng%state = 88172645463393265_i8
    else
      rng%state = seed
    end if
    rng%has_spare = .false.
  end subroutine rng_seed

  function next_u64(rng) result(x)
    type(rng_state), intent(inout) :: rng
    integer(i8) :: x
    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
  end function next_u64

  function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(i8) :: x
    x = next_u64(rng)
    u = real(iand(x, int(z'7FFFFFFFFFFFFFFF', i8)), dp) / real(huge(0_i8), dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function rng_uniform

  function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z, u1, u2, r
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    r = sqrt(-2.0_dp*log(u1))
    z = r*cos(2.0_dp*pi*u2)
    rng%spare = r*sin(2.0_dp*pi*u2)
    rng%has_spare = .true.
  end function rng_normal

  function rng_exponential(rng, mean) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in), optional :: mean
    real(dp) :: x, m
    m = 1.0_dp
    if (present(mean)) m = mean
    x = -m*log(rng_uniform(rng))
  end function rng_exponential

  recursive function rng_gamma(rng, shape, scale) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    real(dp) :: x, sc, d, c, z, u
    sc = 1.0_dp
    if (present(scale)) sc = scale
    if (shape <= 0.0_dp .or. sc < 0.0_dp) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (shape < 1.0_dp) then
      x = rng_gamma(rng, shape + 1.0_dp, sc) * rng_uniform(rng)**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = rng_normal(rng)
        if (1.0_dp + c*z > 0.0_dp) exit
      end do
      u = rng_uniform(rng)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-(1.0_dp+c*z)**3 + log((1.0_dp+c*z)**3))) exit
    end do
    x = sc*d*(1.0_dp+c*z)**3
  end function rng_gamma

end module ecd_rng
