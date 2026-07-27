! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_splines
  use yc_kinds, only : dp
  use yc_linalg, only : solve_linear
  implicit none
  private
  public :: fit_cubic_spline_coefficients, evaluate_cubic_spline

contains

  subroutine fit_cubic_spline_coefficients(x, y, method, b, c, d, ok)
    real(dp), intent(in) :: x(:), y(:)
    character(len=*), intent(in) :: method
    real(dp), allocatable, intent(out) :: b(:), c(:), d(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: amat(:,:), rhs(:), second(:), h(:)
    integer :: n, i
    character(len=:), allocatable :: meth

    n = size(x)
    ok = .false.
    allocate(b(max(0,n-1)), c(max(0,n-1)), d(max(0,n-1)))
    if (n < 3 .or. size(y) /= n) return
    allocate(h(n-1))
    h = x(2:n) - x(1:n-1)
    if (any(h <= 0.0_dp)) return
    allocate(amat(n,n), rhs(n))
    amat = 0.0_dp
    rhs = 0.0_dp
    meth = trim(adjustl(method))

    do i = 2, n - 1
      amat(i,i-1) = h(i-1)
      amat(i,i) = 2.0_dp * (h(i-1) + h(i))
      amat(i,i+1) = h(i)
      rhs(i) = 6.0_dp * ((y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1))
    end do

    if (meth == 'fmm' .and. n >= 4) then
      amat(1,1) = -h(2)
      amat(1,2) = h(1) + h(2)
      amat(1,3) = -h(1)
      amat(n,n-2) = -h(n-1)
      amat(n,n-1) = h(n-2) + h(n-1)
      amat(n,n) = -h(n-2)
    else
      amat(1,1) = 1.0_dp
      amat(n,n) = 1.0_dp
    end if

    call solve_linear(amat, rhs, second, ok)
    if (.not. ok) return
    do i = 1, n - 1
      b(i) = (y(i+1)-y(i))/h(i) - h(i)*(2.0_dp*second(i)+second(i+1))/6.0_dp
      c(i) = second(i) / 2.0_dp
      d(i) = (second(i+1)-second(i)) / (6.0_dp*h(i))
    end do
    ok = .true.
  end subroutine fit_cubic_spline_coefficients

  pure real(dp) function evaluate_cubic_spline(x, y, b, c, d, xout) result(value)
    real(dp), intent(in) :: x(:), y(:), b(:), c(:), d(:), xout
    real(dp) :: dx
    integer :: i, lo, hi, mid, n
    n = size(x)
    if (xout <= x(1)) then
      i = 1
    else if (xout >= x(n)) then
      i = n - 1
    else
      lo = 1
      hi = n
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (x(mid) <= xout) then
          lo = mid
        else
          hi = mid
        end if
      end do
      i = lo
    end if
    dx = xout - x(i)
    value = y(i) + b(i)*dx + c(i)*dx*dx + d(i)*dx*dx*dx
  end function evaluate_cubic_spline

end module yc_splines
