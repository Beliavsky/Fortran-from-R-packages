! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_random
  use multiatsm_kinds, only : dp, pi
  implicit none
  private
  public :: set_random_seed, random_normal, random_integer

contains

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer, allocatable :: put(:)
    integer :: n, i

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(abs(seed) + 104729 * i + 8191 * i * i, huge(1) - 1) + 1
    end do
    call random_seed(put=put)
  end subroutine set_random_seed

  function random_integer(upper) result(value)
    integer, intent(in) :: upper
    integer :: value
    real(dp) :: u

    if (upper <= 1) then
      value = 1
      return
    end if
    call random_number(u)
    value = min(upper, 1 + int(u * real(upper, dp)))
  end function random_integer

  function random_normal() result(value)
    real(dp) :: value
    real(dp) :: u1, u2

    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
  end function random_normal

end module multiatsm_random
