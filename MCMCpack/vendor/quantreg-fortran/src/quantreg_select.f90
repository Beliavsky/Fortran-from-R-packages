! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_select
  use quantreg_kinds, only : dp
  implicit none
  private
  public :: qselect, kuantiles
contains

  real(dp) function qselect(x, prob) result(q)
    real(dp), intent(in) :: x(:), prob
    real(dp), allocatable :: work(:)
    integer :: k, n

    n = size(x)
    if (n == 0) then
      q = 0.0_dp
      return
    end if
    allocate(work(n))
    work = x
    k = max(1, min(n, nint(max(0.0_dp, min(1.0_dp, prob)) * real(n,dp))))
    call select_k(work, k)
    q = work(k)
  end function qselect

  subroutine kuantiles(x, probs, values, qtype)
    real(dp), intent(in) :: x(:), probs(:)
    real(dp), intent(out) :: values(:)
    integer, intent(in), optional :: qtype
    real(dp), allocatable :: s(:)
    real(dp) :: p, a, b, d, g, h
    integer :: i, j, k, n, typ

    n = size(x)
    typ = 7
    if (present(qtype)) typ = qtype
    if (size(values) /= size(probs) .or. n == 0 .or. typ < 1 .or. typ > 9) then
      if (size(values) > 0) values = 0.0_dp
      return
    end if
    allocate(s(n))
    s = x
    call sort_real(s)

    do i = 1, size(probs)
      p = max(0.0_dp, min(1.0_dp, probs(i)))
      select case (typ)
      case (1)
        j = max(1, ceiling(p * real(n,dp)))
        values(i) = s(min(j,n))
      case (2)
        j = max(1, floor(p * real(n,dp)))
        k = min(j + 1, n)
        if (p * real(n,dp) > real(j,dp)) then
          g = 1.0_dp
        else
          g = 0.5_dp
        end if
        values(i) = (1.0_dp-g)*s(min(j,n)) + g*s(k)
      case (3)
        j = max(1, nint(p * real(n,dp)))
        values(i) = s(min(j,n))
      case default
        select case (typ)
        case (4)
          a = 0.0_dp
          b = 1.0_dp
        case (5)
          a = 0.5_dp
          b = 0.5_dp
        case (6)
          a = 0.0_dp
          b = 0.0_dp
        case (7)
          a = 1.0_dp
          b = 1.0_dp
        case (8)
          a = 1.0_dp/3.0_dp
          b = 1.0_dp/3.0_dp
        case (9)
          a = 3.0_dp/8.0_dp
          b = 3.0_dp/8.0_dp
        end select
        d = a + p*(1.0_dp-a-b)
        h = p*real(n,dp) + d
        j = floor(h)
        g = h - real(j,dp)
        j = max(1,min(n,j))
        k = max(1,min(n,j+1))
        values(i) = (1.0_dp-g)*s(j) + g*s(k)
      end select
    end do
  end subroutine kuantiles

  subroutine select_k(a, k)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: k
    integer :: left, right, i, j
    real(dp) :: pivot, tmp

    left = 1
    right = size(a)
    do while (left < right)
      pivot = a((left+right)/2)
      i = left
      j = right
      do
        do while (a(i) < pivot)
          i = i + 1
        end do
        do while (a(j) > pivot)
          j = j - 1
        end do
        if (i > j) exit
        tmp = a(i)
        a(i) = a(j)
        a(j) = tmp
        i = i + 1
        j = j - 1
      end do
      if (k <= j) then
        right = j
      else if (k >= i) then
        left = i
      else
        return
      end if
    end do
  end subroutine select_k

  subroutine sort_real(a)
    real(dp), intent(inout) :: a(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(a)
      key = a(i)
      j = i - 1
      do while (j >= 1)
        if (a(j) <= key) exit
        a(j+1) = a(j)
        j = j - 1
      end do
      a(j+1) = key
    end do
  end subroutine sort_real

end module quantreg_select
