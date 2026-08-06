! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_random
  use rrcov_kinds, only : dp
  implicit none
  private
  public :: seed_random, random_normal, random_subset, random_unit_vector
contains
  subroutine seed_random(seed)
    integer, intent(in), optional :: seed
    integer, allocatable :: values(:)
    integer :: i, n, base
    call random_seed(size=n)
    allocate(values(n))
    base = 104729
    if (present(seed)) base = seed
    do i = 1, n
      values(i) = modulo(base + 7919 * i, huge(1) - 1)
      if (values(i) <= 0) values(i) = i
    end do
    call random_seed(put=values)
  end subroutine seed_random

  function random_normal() result(value)
    real(dp) :: value
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
  end function random_normal

  subroutine random_subset(n, k, index)
    integer, intent(in) :: n, k
    integer, intent(out) :: index(k)
    integer, allocatable :: pool(:)
    real(dp) :: u
    integer :: i, j, tmp
    allocate(pool(n))
    pool = [(i, i=1, n)]
    do i = 1, min(k, n)
      call random_number(u)
      j = i + int(u * real(n - i + 1, dp))
      j = min(n, max(i, j))
      tmp = pool(i)
      pool(i) = pool(j)
      pool(j) = tmp
      index(i) = pool(i)
    end do
  end subroutine random_subset

  subroutine random_unit_vector(x)
    real(dp), intent(out) :: x(:)
    real(dp) :: norm
    integer :: i
    do i = 1, size(x)
      x(i) = random_normal()
    end do
    norm = sqrt(sum(x * x))
    if (norm <= tiny(1.0_dp)) then
      x = 0.0_dp
      x(1) = 1.0_dp
    else
      x = x / norm
    end if
  end subroutine random_unit_vector
end module rrcov_random
