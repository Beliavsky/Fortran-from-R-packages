! SPDX-License-Identifier: GPL-3.0-or-later
module rrcov_sort
  use rrcov_kinds, only : dp
  implicit none
  private
  public :: sort_real, sort_real_with_index, order_smallest, rank_values
contains
  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
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
  end subroutine sort_real

  subroutine sort_real_with_index(x, index)
    real(dp), intent(inout) :: x(:)
    integer, intent(inout) :: index(:)
    integer :: i, j, ikey
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      ikey = index(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j + 1) = x(j)
        index(j + 1) = index(j)
        j = j - 1
      end do
      x(j + 1) = key
      index(j + 1) = ikey
    end do
  end subroutine sort_real_with_index

  subroutine order_smallest(x, k, index)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: k
    integer, intent(out) :: index(k)
    real(dp), allocatable :: work(:)
    integer, allocatable :: order(:)
    integer :: i
    allocate(work(size(x)), order(size(x)))
    work = x
    order = [(i, i=1, size(x))]
    call sort_real_with_index(work, order)
    index = order(1:k)
  end subroutine order_smallest

  subroutine rank_values(x, ranks)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: ranks(size(x))
    real(dp), allocatable :: work(:)
    integer, allocatable :: order(:)
    integer :: i, j, n
    real(dp) :: average_rank
    n = size(x)
    allocate(work(n), order(n))
    work = x
    order = [(i, i=1, n)]
    call sort_real_with_index(work, order)
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (abs(work(j + 1) - work(i)) > epsilon(1.0_dp) * &
          max(1.0_dp, abs(work(j + 1)), abs(work(i)))) exit
        j = j + 1
      end do
      average_rank = 0.5_dp * real(i + j, dp)
      ranks(order(i:j)) = average_rank
      i = j + 1
    end do
  end subroutine rank_values
end module rrcov_sort
