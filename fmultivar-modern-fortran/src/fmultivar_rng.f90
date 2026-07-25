! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_rng
  use fmultivar_kinds, only : dp, i8, pi
  implicit none
  private
  public :: rng_state, seed_rng, uniform_rng, normal_rng, gamma_rng, chi_square_rng

  type :: rng_state
    integer(i8) :: s = 88172645463325252_i8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine seed_rng(state, seed)
    type(rng_state), intent(out) :: state
    integer(i8), intent(in) :: seed
    state%s = merge(seed, 88172645463325252_i8, seed /= 0_i8)
    state%has_spare = .false.
    state%spare = 0.0_dp
  end subroutine seed_rng

  function next_u64(state) result(x)
    type(rng_state), intent(inout) :: state
    integer(i8) :: x
    x = state%s
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    state%s = x
  end function next_u64

  function uniform_rng(state) result(u)
    type(rng_state), intent(inout) :: state
    real(dp) :: u
    integer(i8) :: x
    x = next_u64(state)
    u = real(iand(shiftr(x, 11), int(z'001FFFFFFFFFFFFF', i8)), dp) / 9007199254740992.0_dp
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
  end function uniform_rng

  function normal_rng(state) result(z)
    type(rng_state), intent(inout) :: state
    real(dp) :: z, u1, u2, r
    if (state%has_spare) then
      z = state%spare
      state%has_spare = .false.
      return
    end if
    u1 = uniform_rng(state)
    u2 = uniform_rng(state)
    r = sqrt(-2.0_dp*log(u1))
    z = r*cos(2.0_dp*pi*u2)
    state%spare = r*sin(2.0_dp*pi*u2)
    state%has_spare = .true.
  end function normal_rng

  recursive function gamma_rng(state, shape, scale) result(x)
    type(rng_state), intent(inout) :: state
    real(dp), intent(in) :: shape, scale
    real(dp) :: x, d, c, z, v, u
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = gamma_rng(state, shape + 1.0_dp, scale) * uniform_rng(state)**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      z = normal_rng(state)
      v = 1.0_dp + c*z
      if (v <= 0.0_dp) cycle
      v = v*v*v
      u = uniform_rng(state)
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp - v + log(v))) exit
    end do
    x = scale*d*v
  end function gamma_rng

  function chi_square_rng(state, df) result(x)
    type(rng_state), intent(inout) :: state
    real(dp), intent(in) :: df
    real(dp) :: x
    x = gamma_rng(state, 0.5_dp*df, 2.0_dp)
  end function chi_square_rng
end module fmultivar_rng
