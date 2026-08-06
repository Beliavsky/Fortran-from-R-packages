! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rpeglmen_kinds, only : dp
  implicit none
  private

  public :: safe_exp, vector_norm2, sample_sd, digamma_value
  public :: all_finite, approximately_constant_one

contains

  elemental real(dp) function safe_exp(x) result(value)
    real(dp), intent(in) :: x

    value = exp(max(-700.0_dp, min(700.0_dp, x)))
  end function safe_exp

  pure real(dp) function vector_norm2(x) result(value)
    real(dp), intent(in) :: x(:)

    value = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vector_norm2

  pure real(dp) function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: center

    if (size(x) <= 1) then
      value = 0.0_dp
      return
    end if
    center = sum(x) / real(size(x), dp)
    value = sqrt(max(0.0_dp, sum((x - center)**2) / real(size(x) - 1, dp)))
  end function sample_sd

  pure real(dp) function digamma_value(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: y, inv, inv2

    if (x <= 0.0_dp) then
      value = huge(1.0_dp)
      return
    end if

    y = x
    value = 0.0_dp
    do while (y < 8.0_dp)
      value = value - 1.0_dp / y
      y = y + 1.0_dp
    end do

    inv = 1.0_dp / y
    inv2 = inv * inv
    value = value + log(y) - 0.5_dp * inv &
      - inv2 * (1.0_dp / 12.0_dp &
      - inv2 * (1.0_dp / 120.0_dp &
      - inv2 * (1.0_dp / 252.0_dp &
      - inv2 * (1.0_dp / 240.0_dp &
      - inv2 * (5.0_dp / 660.0_dp)))))
  end function digamma_value

  pure logical function all_finite(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i

    ok = .true.
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) then
        ok = .false.
        return
      end if
    end do
  end function all_finite

  pure logical function approximately_constant_one(x) result(ok)
    real(dp), intent(in) :: x(:)

    if (size(x) == 0) then
      ok = .false.
    else
      ok = maxval(abs(x - 1.0_dp)) <= 100.0_dp * epsilon(1.0_dp)
    end if
  end function approximately_constant_one

end module rpeglmen_math
