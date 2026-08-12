! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_utils
  use adagio_kinds, only : dp
  implicit none
  private
  public :: argsort_real, sort_real, argsort_int, norm2_vec, outer_product

contains

  subroutine argsort_real(x, idx, ascending)
    real(dp), intent(in) :: x(:)
    integer, intent(out) :: idx(size(x))
    logical, intent(in), optional :: ascending
    logical :: asc
    integer :: i, j, key
    asc = .true.; if (present(ascending)) asc = ascending
    do i = 1, size(x); idx(i) = i; end do
    do i = 2, size(x)
       key = idx(i); j = i - 1
       if (asc) then
          do while (j >= 1)
             if (x(idx(j)) <= x(key)) exit
             idx(j+1) = idx(j); j = j - 1
          end do
       else
          do while (j >= 1)
             if (x(idx(j)) >= x(key)) exit
             idx(j+1) = idx(j); j = j - 1
          end do
       end if
       idx(j+1) = key
    end do
  end subroutine argsort_real

  subroutine argsort_int(x, idx, ascending)
    integer, intent(in) :: x(:)
    integer, intent(out) :: idx(size(x))
    logical, intent(in), optional :: ascending
    logical :: asc
    integer :: i, j, key
    asc = .true.; if (present(ascending)) asc = ascending
    do i = 1, size(x); idx(i) = i; end do
    do i = 2, size(x)
       key = idx(i); j = i - 1
       if (asc) then
          do while (j >= 1)
             if (x(idx(j)) <= x(key)) exit
             idx(j+1) = idx(j); j = j - 1
          end do
       else
          do while (j >= 1)
             if (x(idx(j)) >= x(key)) exit
             idx(j+1) = idx(j); j = j - 1
          end do
       end if
       idx(j+1) = key
    end do
  end subroutine argsort_int

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: idx(size(x))
    real(dp) :: tmp(size(x))
    call argsort_real(x, idx)
    tmp = x(idx)
    x = tmp
  end subroutine sort_real

  pure function norm2_vec(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v = sqrt(sum(x*x))
  end function norm2_vec

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x), size(y))
    integer :: j
    do j = 1, size(y)
       a(:,j) = x * y(j)
    end do
  end function outer_product
end module adagio_utils
