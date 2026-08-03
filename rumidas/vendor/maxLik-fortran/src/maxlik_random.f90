! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_random
  use, intrinsic :: iso_fortran_env, only: int64
  use maxlik_kinds, only: dp
  implicit none
  private

  type, public :: maxlik_rng
    integer(int64) :: state = 88172645463393265_int64
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
    procedure :: integer => rng_integer
    procedure :: shuffle => rng_shuffle
  end type maxlik_rng

contains

  subroutine rng_seed(self, seed)
    class(maxlik_rng), intent(inout) :: self
    integer, intent(in) :: seed
    self%state = int(seed, int64)
    if (self%state == 0_int64) self%state = 88172645463393265_int64
    self%state = ieor(self%state, int(z'9E3779B97F4A7C15', int64))
  end subroutine rng_seed

  real(dp) function rng_uniform(self) result(value)
    class(maxlik_rng), intent(inout) :: self
    integer(int64) :: x, positive
    x = self%state
    x = ieor(x, shiftl(x, 13))
    x = ieor(x, shiftr(x, 7))
    x = ieor(x, shiftl(x, 17))
    self%state = x
    positive = iand(x, int(z'7FFFFFFFFFFFFFFF', int64))
    value = (real(positive, dp) + 0.5_dp) / (real(huge(positive), dp) + 1.0_dp)
    value = min(1.0_dp - epsilon(1.0_dp), max(tiny(1.0_dp), value))
  end function rng_uniform

  real(dp) function rng_normal(self) result(value)
    class(maxlik_rng), intent(inout) :: self
    real(dp) :: u1, u2
    u1 = self%uniform()
    u2 = self%uniform()
    value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
  end function rng_normal

  integer function rng_integer(self, low, high) result(value)
    class(maxlik_rng), intent(inout) :: self
    integer, intent(in) :: low, high
    if (high <= low) then
      value = low
    else
      value = low + int(self%uniform() * real(high - low + 1, dp))
      value = min(high, value)
    end if
  end function rng_integer

  subroutine rng_shuffle(self, values)
    class(maxlik_rng), intent(inout) :: self
    integer, intent(inout) :: values(:)
    integer :: i, j, temp
    do i = size(values), 2, -1
      j = self%integer(1, i)
      temp = values(i)
      values(i) = values(j)
      values(j) = temp
    end do
  end subroutine rng_shuffle

end module maxlik_random
