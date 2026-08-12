! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_rng
  use iso_fortran_env, only : int64
  use adagio_kinds, only : dp, pi
  implicit none
  private
  public :: rng_state

  type :: rng_state
     integer(int64) :: state = 88172645463393265_int64
     logical :: have_spare = .false.
     real(dp) :: spare = 0.0_dp
   contains
     procedure :: seed => rng_seed
     procedure :: uniform => rng_uniform
     procedure :: normal => rng_normal
     procedure :: permutation => rng_permutation
  end type rng_state

contains

  subroutine rng_seed(self, seed)
    class(rng_state), intent(inout) :: self
    integer(int64), intent(in) :: seed
    self%state = seed
    if (self%state == 0_int64) self%state = 88172645463393265_int64
    self%have_spare = .false.
  end subroutine rng_seed

  function rng_uniform(self) result(u)
    class(rng_state), intent(inout) :: self
    real(dp) :: u
    integer(int64) :: x
    x = self%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    self%state = x
    u = real(iand(x, int(z'001FFFFFFFFFFFFF', int64)), dp) / real(2_int64**53, dp)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
  end function rng_uniform

  function rng_normal(self) result(z)
    class(rng_state), intent(inout) :: self
    real(dp) :: z, u1, u2, r
    if (self%have_spare) then
       z = self%spare
       self%have_spare = .false.
       return
    end if
    u1 = self%uniform()
    u2 = self%uniform()
    r = sqrt(-2.0_dp * log(u1))
    z = r * cos(2.0_dp*pi*u2)
    self%spare = r * sin(2.0_dp*pi*u2)
    self%have_spare = .true.
  end function rng_normal

  subroutine rng_permutation(self, p)
    class(rng_state), intent(inout) :: self
    integer, intent(out) :: p(:)
    integer :: i, j, t
    do i = 1, size(p)
       p(i) = i
    end do
    do i = size(p), 2, -1
       j = 1 + int(self%uniform() * real(i, dp))
       if (j > i) j = i
       t = p(i); p(i) = p(j); p(j) = t
    end do
  end subroutine rng_permutation

end module adagio_rng
