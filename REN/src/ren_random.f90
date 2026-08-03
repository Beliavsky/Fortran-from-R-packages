! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_random
  use ren_kinds, only : dp
  implicit none
  private
  public :: initialize_random_seed, sample_without_replacement, random_choice
contains
  subroutine initialize_random_seed(seed)
    integer, intent(in) :: seed
    integer, allocatable :: put(:)
    integer :: i, n
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729 * i + 37 * i * i, huge(1) - 1)
      if (put(i) == 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine initialize_random_seed

  subroutine sample_without_replacement(n, k, index)
    integer, intent(in) :: n, k
    integer, intent(out) :: index(k)
    integer, allocatable :: pool(:)
    real(dp) :: u
    integer :: i, j, tmp
    allocate(pool(n))
    pool = [(i, i=1,n)]
    do i = 1, k
      call random_number(u)
      j = i + int(u * real(n - i + 1, dp))
      if (j > n) j = n
      tmp = pool(i)
      pool(i) = pool(j)
      pool(j) = tmp
      index(i) = pool(i)
    end do
  end subroutine sample_without_replacement

  integer function random_choice(index) result(value)
    integer, intent(in) :: index(:)
    real(dp) :: u
    integer :: k
    call random_number(u)
    k = 1 + int(u * real(size(index), dp))
    if (k > size(index)) k = size(index)
    value = index(k)
  end function random_choice
end module ren_random
