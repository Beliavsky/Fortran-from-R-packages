! SPDX-License-Identifier: GPL-3.0-only
module ao_random
  use ao_kinds, only : dp
  use ao_types, only : ao_block
  implicit none
  private
  public :: ao_seed, generate_random_partition
contains
  subroutine ao_seed(seed)
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
  end subroutine ao_seed

  subroutine shuffle_integer(x)
    integer, intent(inout) :: x(:)
    integer :: i, j, tmp
    real(dp) :: u
    do i = size(x), 2, -1
      call random_number(u)
      j = 1 + int(u * real(i, dp))
      if (j > i) j = i
      tmp = x(i); x(i) = x(j); x(j) = tmp
    end do
  end subroutine shuffle_integer

  subroutine sort_integer(x)
    integer, intent(inout) :: x(:)
    integer :: i, j, key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        j = j - 1
      end do
      x(j + 1) = key
    end do
  end subroutine sort_integer

  subroutine generate_random_partition(npar, p, minimum_blocks, blocks)
    integer, intent(in) :: npar, minimum_blocks
    real(dp), intent(in) :: p
    type(ao_block), allocatable, intent(out) :: blocks(:)
    integer, allocatable :: x(:), start_flag(:), zero_pos(:), starts(:)
    integer :: i, j, b, s, e, nz, need
    real(dp) :: u

    if (npar <= 0) then
      allocate(blocks(0)); return
    end if
    if (minimum_blocks >= npar) then
      allocate(blocks(npar))
      do i = 1, npar
        allocate(blocks(i)%index(1)); blocks(i)%index = i
      end do
      return
    end if

    allocate(x(npar), start_flag(npar))
    x = [(i, i=1,npar)]
    call shuffle_integer(x)
    start_flag = 0
    start_flag(1) = 1
    do i = 2, npar
      call random_number(u)
      if (u < p) start_flag(i) = 1
    end do

    b = count(start_flag == 1)
    if (b < minimum_blocks) then
      nz = count(start_flag == 0)
      allocate(zero_pos(nz))
      j = 0
      do i = 1, npar
        if (start_flag(i) == 0) then
          j = j + 1; zero_pos(j) = i
        end if
      end do
      call shuffle_integer(zero_pos)
      need = minimum_blocks - b
      do i = 1, need
        start_flag(zero_pos(i)) = 1
      end do
    end if

    b = count(start_flag == 1)
    allocate(starts(b), blocks(b))
    j = 0
    do i = 1, npar
      if (start_flag(i) == 1) then
        j = j + 1; starts(j) = i
      end if
    end do
    do j = 1, b
      s = starts(j)
      if (j < b) then
        e = starts(j + 1) - 1
      else
        e = npar
      end if
      allocate(blocks(j)%index(e - s + 1))
      blocks(j)%index = x(s:e)
      call sort_integer(blocks(j)%index)
    end do
  end subroutine generate_random_partition
end module ao_random
