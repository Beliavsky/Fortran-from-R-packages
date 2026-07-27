! SPDX-License-Identifier: MIT
! Copyright (c) 2023 Bernardo Reckziegel
module epo_linalg
  use epo_kinds, only : dp
  implicit none
  private

  public :: solve_spd
  public :: quadratic_form

contains

  subroutine solve_spd(a, b, x, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: l(:,:), y(:)
    real(dp) :: value, tol, max_diag
    integer :: i, j, k, n

    n = size(b)
    ok = .false.
    x = 0.0_dp

    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return

    allocate(l(n,n), y(n))
    l = 0.0_dp

    max_diag = 0.0_dp
    do i = 1, n
      max_diag = max(max_diag, abs(a(i,i)))
    end do
    tol = 100.0_dp * epsilon(1.0_dp) * real(max(1,n), dp) * max(1.0_dp, max_diag)

    do i = 1, n
      do j = 1, i
        value = a(i,j)
        do k = 1, j - 1
          value = value - l(i,k) * l(j,k)
        end do

        if (i == j) then
          if (value <= tol) return
          l(i,j) = sqrt(value)
        else
          l(i,j) = value / l(j,j)
        end if
      end do
    end do

    do i = 1, n
      value = b(i)
      do k = 1, i - 1
        value = value - l(i,k) * y(k)
      end do
      y(i) = value / l(i,i)
    end do

    do i = n, 1, -1
      value = y(i)
      do k = i + 1, n
        value = value - l(k,i) * x(k)
      end do
      x(i) = value / l(i,i)
    end do

    ok = .true.
  end subroutine solve_spd

  pure function quadratic_form(x, a) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: value

    value = dot_product(x, matmul(a, x))
  end function quadratic_form

end module epo_linalg
