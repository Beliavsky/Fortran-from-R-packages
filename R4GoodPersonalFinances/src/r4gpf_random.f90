! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_random
  use r4gpf_kinds, only: dp, i8, pi
  implicit none
  private
  public :: rng_state, seed_rng, random_uniform, random_normal, fill_normal

  type :: rng_state
    integer(i8) :: state = 88172645463393265_i8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state
contains

  subroutine seed_rng(rng, seed)
    type(rng_state), intent(inout) :: rng
    integer(i8), intent(in) :: seed
    rng%state = seed
    if (rng%state == 0_i8) rng%state = 88172645463393265_i8
    rng%has_spare = .false.
    rng%spare = 0.0_dp
  end subroutine seed_rng

  function next_u64(rng) result(x)
    type(rng_state), intent(inout) :: rng
    integer(i8) :: x
    x = rng%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    rng%state = x
  end function next_u64

  function random_uniform(rng) result(u)
    type(rng_state), intent(inout) :: rng
    real(dp) :: u
    integer(i8) :: x
    x = next_u64(rng)
    u = real(iand(x, int(z'001FFFFFFFFFFFFF', i8)), dp) / real(int(z'0020000000000000', i8), dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function random_uniform

  function random_normal(rng) result(z)
    type(rng_state), intent(inout) :: rng
    real(dp) :: z, u1, u2, r
    if (rng%has_spare) then
      z = rng%spare
      rng%has_spare = .false.
      return
    end if
    u1 = random_uniform(rng)
    u2 = random_uniform(rng)
    r = sqrt(-2.0_dp * log(u1))
    z = r * cos(2.0_dp * pi * u2)
    rng%spare = r * sin(2.0_dp * pi * u2)
    rng%has_spare = .true.
  end function random_normal

  subroutine fill_normal(rng, x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(out) :: x(:)
    integer :: i
    do i = 1, size(x)
      x(i) = random_normal(rng)
    end do
  end subroutine fill_normal

end module r4gpf_random
