! SPDX-License-Identifier: GPL-2.0-or-later
module mcga_random
  use mcga_kinds, only : dp
  implicit none
  private
  public :: set_random_seed, runif_scalar, randint
contains
  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)

    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729 * (i - 1), huge(1) - 1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine set_random_seed

  function runif_scalar(a, b) result(x)
    real(dp), intent(in) :: a, b
    real(dp) :: x, u
    call random_number(u)
    x = a + (b - a) * u
  end function runif_scalar

  function randint(lo, hi) result(k)
    integer, intent(in) :: lo, hi
    integer :: k
    real(dp) :: u
    if (hi < lo) error stop "randint: invalid range"
    call random_number(u)
    k = lo + int(u * real(hi - lo + 1, dp))
    if (k > hi) k = hi
  end function randint
end module mcga_random
