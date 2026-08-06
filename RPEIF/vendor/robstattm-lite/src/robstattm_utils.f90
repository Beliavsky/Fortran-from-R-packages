! SPDX-License-Identifier: GPL-3.0-or-later
!
! Minimal utility subset required by robstattm_psi.f90.  The interfaces and
! behavior follow the completed RobStatTM modern Fortran translation.
module robstattm_utils
  use robstattm_kinds, only : dp
  implicit none
  private
  public :: lower_string, mean_value, median_absolute
contains
  pure function lower_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, code

    out = text
    do i = 1, len(text)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
    end do
  end function lower_string

  pure function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value

    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  function median_absolute(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    real(dp), allocatable :: work(:)
    integer :: n

    n = size(x)
    if (n == 0) then
      value = 0.0_dp
      return
    end if
    allocate(work(n))
    work = abs(x)
    call sort_real(work)
    if (mod(n, 2) == 1) then
      value = work((n + 1) / 2)
    else
      value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
    end if
  end function median_absolute

  recursive subroutine quicksort(a, left, right)
    real(dp), intent(inout) :: a(:)
    integer, intent(in) :: left, right
    integer :: i, j
    real(dp) :: pivot, tmp

    if (left >= right) return
    i = left
    j = right
    pivot = a((left + right) / 2)
    do
      do while (a(i) < pivot)
        i = i + 1
      end do
      do while (a(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tmp = a(i)
        a(i) = a(j)
        a(j) = tmp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (left < j) call quicksort(a, left, j)
    if (i < right) call quicksort(a, i, right)
  end subroutine quicksort

  subroutine sort_real(a)
    real(dp), intent(inout) :: a(:)
    if (size(a) > 1) call quicksort(a, 1, size(a))
  end subroutine sort_real
end module robstattm_utils
