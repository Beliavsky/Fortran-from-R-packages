! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_rng
  use fextremes_kinds, only: dp, pi
  implicit none
  private
  public :: rng_state, seed_rng, uniform_rng, normal_rng, exponential_rng

  type :: rng_state
    integer(kind=8) :: state = 88172645463325252_8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains
  subroutine seed_rng(rng, seed)
    type(rng_state), intent(out) :: rng
    integer, intent(in) :: seed
    rng%state = int(max(1, seed), kind=8)
    rng%state = ieor(rng%state, 4101842887655102017_8)
    rng%has_spare = .false.
  end subroutine seed_rng

  real(dp) function uniform_rng(rng) result(u)
    type(rng_state), intent(inout) :: rng
    integer(kind=8) :: x
    x = rng%state
    x = ieor(x, shiftr(x, 12))
    x = ieor(x, shiftl(x, 25))
    x = ieor(x, shiftr(x, 27))
    rng%state = x
    x = x * 2685821657736338717_8
    u = real(iand(shiftr(x, 11), int(z'001FFFFFFFFFFFFF', kind=8)), dp) / 9007199254740992.0_dp
    u = min(max(u, epsilon(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
  end function uniform_rng

  real(dp) function exponential_rng(rng) result(x)
    type(rng_state), intent(inout) :: rng
    x = -log(uniform_rng(rng))
  end function exponential_rng

  real(dp) function normal_rng(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u1, u2, r
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = uniform_rng(rng)
    u2 = uniform_rng(rng)
    r = sqrt(-2.0_dp * log(u1))
    z = r * cos(2.0_dp * pi * u2)
    rng%spare = r * sin(2.0_dp * pi * u2)
    rng%has_spare = .true.
  end function normal_rng
end module fextremes_rng
