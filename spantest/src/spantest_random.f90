! SPDX-License-Identifier: GPL-3.0-only
module spantest_random
  use iso_fortran_env, only : int64
  use spantest_kinds, only : dp, pi_dp
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_student_t

  type :: rng_state
    integer(int64) :: state = 88172645463325252_int64
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

contains

  subroutine rng_seed(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: seed
    integer(int64) :: s
    s = int(seed, int64)
    if (s == 0_int64) s = 1_int64
    rng%state = ieor(s, int(z'9E3779B97F4A7C15', int64))
    if (rng%state == 0_int64) rng%state = 88172645463325252_int64
    rng%has_spare = .false.
  end subroutine rng_seed

  function rng_uint64(rng) result(x)
    type(rng_state), intent(inout) :: rng
    integer(int64) :: x
    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
  end function rng_uint64

  real(dp) function rng_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    integer(int64) :: x
    x = rng_uint64(rng)
    u = real(iand(shiftr(x,11), int(z'001FFFFFFFFFFFFF',int64)), dp) / &
        9007199254740992.0_dp
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function rng_uniform

  real(dp) function rng_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1, u2, r
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = rng_uniform(rng)
    u2 = rng_uniform(rng)
    r = sqrt(-2.0_dp*log(u1))
    z = r*cos(2.0_dp*pi_dp*u2)
    rng%spare = r*sin(2.0_dp*pi_dp*u2)
    rng%has_spare = .true.
  end function rng_normal

  recursive real(dp) function rng_gamma(rng, shape) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: shape
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp) then
      x = 0.0_dp
    else if (shape < 1.0_dp) then
      x = rng_gamma(rng, shape+1.0_dp) * rng_uniform(rng)**(1.0_dp/shape)
    else
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
        do
          z = rng_normal(rng)
          v = 1.0_dp + c*z
          if (v > 0.0_dp) exit
        end do
        v = v*v*v
        u = rng_uniform(rng)
        if (u < 1.0_dp - 0.0331_dp*z**4) exit
        if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
      end do
      x = d*v
    end if
  end function rng_gamma

  real(dp) function rng_student_t(rng, df) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: df
    real(dp) :: chi2
    chi2 = 2.0_dp*rng_gamma(rng, 0.5_dp*df)
    x = rng_normal(rng)/sqrt(chi2/df)
  end function rng_student_t

end module spantest_random
