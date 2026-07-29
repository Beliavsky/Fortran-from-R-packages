! SPDX-License-Identifier: GPL-2.0-or-later
module test_support
  use apt_kinds, only : dp
  implicit none
  private
  public :: assert_close, assert_true, generate_prices
contains
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual-expected) > tolerance * max(1.0_dp, abs(expected))) then
      write(*,'(a)') 'FAIL: '//trim(message)
      write(*,'(a,es24.15)') '  actual:   ', actual
      write(*,'(a,es24.15)') '  expected: ', expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine generate_prices(x, y)
    real(dp), intent(out) :: x(:), y(:)
    real(dp), allocatable :: z(:)
    real(dp) :: dx, dzlag, phi
    integer :: t, n
    n = size(x)
    if (size(y) /= n) error stop 'generate_prices: dimensions differ'
    allocate(z(n))
    x(1) = 10.0_dp
    z(1) = 0.8_dp
    do t = 2, n
      dx = 0.12_dp + 0.22_dp*sin(0.37_dp*real(t,dp)) + &
        0.08_dp*cos(0.11_dp*real(t,dp))
      x(t) = x(t-1) + dx
      if (t == 2) then
        dzlag = 0.0_dp
      else
        dzlag = z(t-1)-z(t-2)
      end if
      if (z(t-1) >= 0.0_dp) then
        phi = -0.18_dp
      else
        phi = -0.42_dp
      end if
      z(t) = z(t-1) + phi*z(t-1) + 0.22_dp*dzlag + &
        0.07_dp*sin(0.73_dp*real(t,dp))
    end do
    do t = 1, n
      y(t) = 2.5_dp + 1.35_dp*x(t) + z(t) + &
        0.03_dp*cos(0.21_dp*real(t,dp))
    end do
  end subroutine generate_prices
end module test_support
