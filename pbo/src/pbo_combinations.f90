! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo_combinations
  use iso_fortran_env, only : int64
  implicit none
  private
  public :: binomial_coefficient, generate_combinations
contains
  pure function binomial_coefficient(n, k) result(value)
    integer, intent(in) :: n, k
    integer(int64) :: value
    integer :: i, kk

    if (k < 0 .or. k > n) then
      value = 0_int64
      return
    end if
    kk = min(k, n - k)
    value = 1_int64
    do i = 1, kk
      value = value * int(n - kk + i, int64) / int(i, int64)
    end do
  end function binomial_coefficient

  subroutine generate_combinations(n, k, combinations, success, message)
    integer, intent(in) :: n, k
    integer, allocatable, intent(out) :: combinations(:,:)
    logical, intent(out) :: success
    character(len=:), allocatable, intent(out) :: message
    integer(int64) :: nc64
    integer :: nc, col, i, j
    integer, allocatable :: current(:)

    success = .false.
    message = ''
    if (n < 1 .or. k < 0 .or. k > n) then
      message = 'invalid combination dimensions'
      allocate(combinations(0,0))
      return
    end if
    nc64 = binomial_coefficient(n, k)
    if (nc64 > int(huge(1), int64)) then
      message = 'number of combinations exceeds default integer capacity'
      allocate(combinations(0,0))
      return
    end if
    nc = int(nc64)
    allocate(combinations(k,nc), current(k))
    if (k == 0) then
      success = .true.
      return
    end if
    do i = 1, k
      current(i) = i
    end do
    do col = 1, nc
      combinations(:,col) = current
      i = k
      do while (i >= 1)
        if (current(i) < n - k + i) exit
        i = i - 1
      end do
      if (i < 1) exit
      current(i) = current(i) + 1
      do j = i + 1, k
        current(j) = current(j - 1) + 1
      end do
    end do
    success = .true.
  end subroutine generate_combinations
end module pbo_combinations
