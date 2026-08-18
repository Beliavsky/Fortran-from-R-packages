! SPDX-License-Identifier: MIT
module chyper_rng
  use chyper_kinds, only : dp
  use chyper_distribution, only : qchyper
  implicit none
  private
  public :: chyper_seed, rchyper, rchyper_one
contains
  subroutine chyper_seed(seed)
    integer, intent(in) :: seed
    integer, allocatable :: put(:)
    integer :: n, i
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729 * (i - 1), huge(1) - 1) + 1
    end do
    call random_seed(put=put)
  end subroutine chyper_seed

  integer function rchyper_one(s, n, m) result(x)
    integer, intent(in) :: s
    integer, intent(in) :: n(:), m(:)
    real(dp) :: u
    call random_number(u)
    x = qchyper(u, s, n, m)
  end function rchyper_one

  subroutine rchyper(size_out, s, n, m, x)
    integer, intent(in) :: size_out, s
    integer, intent(in) :: n(:), m(:)
    integer, intent(out) :: x(size_out)
    integer :: i
    do i = 1, size_out
      x(i) = rchyper_one(s, n, m)
    end do
  end subroutine rchyper
end module chyper_rng
