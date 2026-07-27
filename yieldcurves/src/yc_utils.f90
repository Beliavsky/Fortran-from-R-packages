! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_utils
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use yc_kinds, only : dp
  implicit none
  private
  public :: valid_positive_vector, valid_finite_vector, sort_pairs, linear_interpolate
  public :: lower_string, quiet_nan, linspace, payment_times

contains

  pure logical function valid_positive_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    ok = size(x) > 0
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i)) .or. x(i) <= 0.0_dp) then
        ok = .false.
        return
      end if
    end do
  end function valid_positive_vector

  pure logical function valid_finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    ok = size(x) > 0
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) then
        ok = .false.
        return
      end if
    end do
  end function valid_finite_vector

  subroutine sort_pairs(x, y, w)
    real(dp), intent(inout) :: x(:), y(:)
    real(dp), intent(inout), optional :: w(:)
    real(dp) :: tx, ty, tw
    integer :: i, j
    do i = 2, size(x)
      tx = x(i)
      ty = y(i)
      if (present(w)) tw = w(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= tx) exit
        x(j+1) = x(j)
        y(j+1) = y(j)
        if (present(w)) w(j+1) = w(j)
        j = j - 1
      end do
      x(j+1) = tx
      y(j+1) = ty
      if (present(w)) w(j+1) = tw
    end do
  end subroutine sort_pairs

  pure real(dp) function linear_interpolate(x, y, xout) result(value)
    real(dp), intent(in) :: x(:), y(:), xout
    integer :: lo, hi, mid
    if (xout <= x(1)) then
      value = y(1)
      return
    end if
    if (xout >= x(size(x))) then
      value = y(size(y))
      return
    end if
    lo = 1
    hi = size(x)
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x(mid) <= xout) then
        lo = mid
      else
        hi = mid
      end if
    end do
    value = y(lo) + (y(hi) - y(lo)) * (xout - x(lo)) / (x(hi) - x(lo))
  end function linear_interpolate

  pure function lower_string(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, c
    lower = text
    do i = 1, len(text)
      c = iachar(text(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) lower(i:i) = achar(c + 32)
    end do
  end function lower_string

  real(dp) function quiet_nan() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function quiet_nan

  function linspace(a, b, n) result(x)
    real(dp), intent(in) :: a, b
    integer, intent(in) :: n
    real(dp), allocatable :: x(:)
    integer :: i
    allocate(x(n))
    if (n == 1) then
      x(1) = a
    else
      do i = 1, n
        x(i) = a + (b-a) * real(i-1,dp) / real(n-1,dp)
      end do
    end if
  end function linspace

  function payment_times(maturity, frequency, ok) result(times)
    real(dp), intent(in) :: maturity
    integer, intent(in) :: frequency
    logical, intent(out) :: ok
    real(dp), allocatable :: times(:)
    integer :: n, i
    ok = .false.
    allocate(times(0))
    if (maturity <= 0.0_dp) return
    if (frequency /= 1 .and. frequency /= 2) return
    n = nint(maturity * real(frequency,dp))
    if (n < 1) return
    if (abs(real(n,dp)/real(frequency,dp)-maturity) > 1.0e-8_dp) return
    deallocate(times)
    allocate(times(n))
    do i = 1, n
      times(i) = real(i,dp) / real(frequency,dp)
    end do
    times(n) = maturity
    ok = .true.
  end function payment_times

end module yc_utils
