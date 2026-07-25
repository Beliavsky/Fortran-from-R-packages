! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear_rng
  use chaos_kinds, only : dp
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer(int64), parameter :: modulus = 2147483647_int64
  integer(int64), parameter :: multiplier = 16807_int64
  integer(int64), parameter :: quotient = 127773_int64
  integer(int64), parameter :: remainder = 2836_int64

  type, public :: rng_state
    integer(int64) :: seed = 123456789_int64
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

  public :: rng_seed, rng_uniform, rng_normal, fill_uniform, fill_normal
contains
  subroutine rng_seed(state, seed)
    type(rng_state), intent(inout) :: state
    integer(int64), intent(in) :: seed
    state%seed = modulo(abs(seed), modulus - 1_int64) + 1_int64
    state%has_spare = .false.
    state%spare = 0.0_dp
  end subroutine rng_seed

  real(dp) function rng_uniform(state) result(u)
    type(rng_state), intent(inout) :: state
    integer(int64) :: high, low, test
    high = state%seed / quotient
    low = modulo(state%seed, quotient)
    test = multiplier * low - remainder * high
    if (test > 0_int64) then
      state%seed = test
    else
      state%seed = test + modulus
    end if
    u = real(state%seed, dp) / real(modulus, dp)
  end function rng_uniform

  real(dp) function rng_normal(state) result(z)
    type(rng_state), intent(inout) :: state
    real(dp) :: u1, u2, radius
    real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp
    if (state%has_spare) then
      z = state%spare
      state%has_spare = .false.
      return
    end if
    u1 = max(rng_uniform(state), tiny(1.0_dp))
    u2 = rng_uniform(state)
    radius = sqrt(-2.0_dp * log(u1))
    z = radius * cos(two_pi * u2)
    state%spare = radius * sin(two_pi * u2)
    state%has_spare = .true.
  end function rng_normal

  subroutine fill_uniform(state, x, lower, upper)
    type(rng_state), intent(inout) :: state
    real(dp), intent(out) :: x(..)
    real(dp), intent(in), optional :: lower, upper
    real(dp) :: lo, hi
    integer :: i, j
    lo = 0.0_dp
    hi = 1.0_dp
    if (present(lower)) lo = lower
    if (present(upper)) hi = upper
    select rank (x)
    rank (1)
      do i = 1, size(x)
        x(i) = lo + (hi - lo) * rng_uniform(state)
      end do
    rank (2)
      do j = 1, size(x, 2)
        do i = 1, size(x, 1)
          x(i, j) = lo + (hi - lo) * rng_uniform(state)
        end do
      end do
    rank default
      error stop "fill_uniform: unsupported rank"
    end select
  end subroutine fill_uniform

  subroutine fill_normal(state, x)
    type(rng_state), intent(inout) :: state
    real(dp), intent(out) :: x(..)
    integer :: i, j
    select rank (x)
    rank (1)
      do i = 1, size(x)
        x(i) = rng_normal(state)
      end do
    rank (2)
      do j = 1, size(x, 2)
        do i = 1, size(x, 1)
          x(i, j) = rng_normal(state)
        end do
      end do
    rank default
      error stop "fill_normal: unsupported rank"
    end select
  end subroutine fill_normal
end module fnonlinear_rng
