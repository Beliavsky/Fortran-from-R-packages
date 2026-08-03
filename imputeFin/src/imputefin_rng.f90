! SPDX-License-Identifier: GPL-3.0-only
module imputefin_rng
  use imputefin_kinds, only : dp, pi
  implicit none
  private
  public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq, rng_student_t

  type :: rng_state
    integer(kind=8) :: s = 5489_8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine rng_seed(state, seed)
    type(rng_state), intent(out) :: state
    integer(kind=8), intent(in) :: seed
    state%s = merge(seed, 5489_8, seed /= 0_8)
    state%has_spare = .false.
  end subroutine rng_seed

  function next_u64(state) result(x)
    type(rng_state), intent(inout) :: state
    integer(kind=8) :: x
    x = state%s
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    state%s = x
  end function next_u64

  function rng_uniform(state) result(u)
    type(rng_state), intent(inout) :: state
    real(dp) :: u
    integer(kind=8) :: x
    x = next_u64(state)
    u = real(iand(shiftr(x, 11), int(z'001FFFFFFFFFFFFF', kind=8)), dp) / 9007199254740992.0_dp
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
  end function rng_uniform

  function rng_normal(state) result(z)
    type(rng_state), intent(inout) :: state
    real(dp) :: z, r, theta
    if (state%has_spare) then
      z = state%spare
      state%has_spare = .false.
    else
      r = sqrt(-2.0_dp * log(rng_uniform(state)))
      theta = 2.0_dp * pi * rng_uniform(state)
      z = r * cos(theta)
      state%spare = r * sin(theta)
      state%has_spare = .true.
    end if
  end function rng_normal

  recursive function rng_gamma(state, shape, rate) result(x)
    type(rng_state), intent(inout) :: state
    real(dp), intent(in) :: shape, rate
    real(dp) :: x, d, c, z, v, u
    if (shape <= 0.0_dp .or. rate <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = rng_gamma(state, shape + 1.0_dp, rate) * rng_uniform(state)**(1.0_dp / shape)
      return
    end if
    d = shape - 1.0_dp / 3.0_dp
    c = 1.0_dp / sqrt(9.0_dp * d)
    do
      z = rng_normal(state)
      v = 1.0_dp + c * z
      if (v <= 0.0_dp) cycle
      v = v**3
      u = rng_uniform(state)
      if (u < 1.0_dp - 0.0331_dp * z**4) exit
      if (log(u) < 0.5_dp * z*z + d * (1.0_dp - v + log(v))) exit
    end do
    x = d * v / rate
  end function rng_gamma

  function rng_chisq(state, df) result(x)
    type(rng_state), intent(inout) :: state
    real(dp), intent(in) :: df
    real(dp) :: x
    x = rng_gamma(state, 0.5_dp * df, 0.5_dp)
  end function rng_chisq

  function rng_student_t(state, df) result(x)
    type(rng_state), intent(inout) :: state
    real(dp), intent(in) :: df
    real(dp) :: x
    x = rng_normal(state) / sqrt(rng_chisq(state, df) / df)
  end function rng_student_t
end module imputefin_rng
