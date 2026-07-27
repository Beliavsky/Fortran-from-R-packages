! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask_rng
  use iso_fortran_env, only: int64
  use bidask_kinds, only: dp
  implicit none
  private

  type, public :: rng_state
    integer(int64) :: state = 88172645463393265_int64
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
    procedure :: bernoulli => rng_bernoulli
  end type rng_state

contains

  subroutine rng_seed(self, seed)
    class(rng_state), intent(inout) :: self
    integer(int64), intent(in) :: seed
    if (seed == 0_int64) then
      self%state = 88172645463393265_int64
    else
      self%state = seed
    end if
    self%has_spare = .false.
    self%spare = 0.0_dp
  end subroutine rng_seed

  real(dp) function rng_uniform(self) result(u)
    class(rng_state), intent(inout) :: self
    integer(int64) :: x, bits
    x = self%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    self%state = x
    bits = iand(x, int(z'001FFFFFFFFFFFFF', int64))
    u = (real(bits, dp) + 0.5_dp) / 9007199254740992.0_dp
  end function rng_uniform

  real(dp) function rng_normal(self) result(z)
    class(rng_state), intent(inout) :: self
    real(dp) :: u1, u2, radius
    real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

    if (self%has_spare) then
      z = self%spare
      self%has_spare = .false.
      return
    end if
    u1 = max(self%uniform(), tiny(1.0_dp))
    u2 = self%uniform()
    radius = sqrt(-2.0_dp * log(u1))
    z = radius * cos(twopi * u2)
    self%spare = radius * sin(twopi * u2)
    self%has_spare = .true.
  end function rng_normal

  logical function rng_bernoulli(self, probability) result(draw)
    class(rng_state), intent(inout) :: self
    real(dp), intent(in) :: probability
    if (probability <= 0.0_dp) then
      draw = .false.
    else if (probability >= 1.0_dp) then
      draw = .true.
    else
      draw = self%uniform() < probability
    end if
  end function rng_bernoulli

end module bidask_rng
