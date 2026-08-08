! SPDX-License-Identifier: GPL-3.0-only
module rgenoud_random
  use rgenoud_kinds, only : dp
  implicit none
  private
  public :: seed_rng, randu, randi, coin_flip
contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729 * i + 37 * i * i, huge(1) - 1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  real(dp) function randu(lo, hi) result(x)
    real(dp), intent(in), optional :: lo, hi
    real(dp) :: u, a, b
    call random_number(u)
    a = 0.0_dp
    b = 1.0_dp
    if (present(lo)) a = lo
    if (present(hi)) b = hi
    x = a + (b - a) * u
  end function randu

  integer function randi(lo, hi) result(k)
    integer, intent(in) :: lo, hi
    real(dp) :: u
    if (hi <= lo) then
      k = lo
      return
    end if
    call random_number(u)
    k = lo + int(u * real(hi - lo + 1, dp))
    if (k > hi) k = hi
  end function randi

  logical function coin_flip() result(head)
    head = (randi(0, 1) == 1)
  end function coin_flip
end module rgenoud_random
