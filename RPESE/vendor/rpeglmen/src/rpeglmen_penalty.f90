! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen_penalty
  use rpeglmen_kinds, only : dp
  implicit none
  private

  public :: prox_l1, prox_en, regularizer_en

contains

  elemental real(dp) function soft_threshold_scalar(x, threshold) result(value)
    real(dp), intent(in) :: x, threshold

    value = sign(max(abs(x) - threshold, 0.0_dp), x)
  end function soft_threshold_scalar

  pure function prox_l1(x, threshold, has_intercept) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: threshold
    logical, intent(in), optional :: has_intercept
    real(dp) :: value(size(x))
    logical :: keep_first
    integer :: j

    keep_first = .false.
    if (present(has_intercept)) keep_first = has_intercept

    do j = 1, size(x)
      value(j) = soft_threshold_scalar(x(j), threshold)
    end do
    if (keep_first .and. size(x) > 0) value(1) = x(1)
  end function prox_l1

  pure function prox_en(x, step, lambda, alpha, penalize_intercept, source_proximal) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: step, lambda, alpha
    logical, intent(in) :: penalize_intercept, source_proximal
    real(dp) :: value(size(x))
    real(dp) :: denominator, threshold
    integer :: j, first_penalized

    value = x
    first_penalized = 1
    if (.not. penalize_intercept .and. size(x) > 0) first_penalized = 2

    denominator = 1.0_dp + step * lambda * (1.0_dp - alpha)
    if (source_proximal) then
      threshold = step * lambda * alpha / denominator
      do j = first_penalized, size(x)
        value(j) = soft_threshold_scalar(x(j), threshold)
      end do
    else
      threshold = step * lambda * alpha
      do j = first_penalized, size(x)
        value(j) = soft_threshold_scalar(x(j), threshold) / denominator
      end do
    end if
  end function prox_en

  pure real(dp) function regularizer_en(x, alpha, penalize_intercept) result(value)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: alpha
    logical, intent(in) :: penalize_intercept
    integer :: first_penalized

    first_penalized = 1
    if (.not. penalize_intercept .and. size(x) > 0) first_penalized = 2

    if (first_penalized > size(x)) then
      value = 0.0_dp
    else
      value = alpha * sum(abs(x(first_penalized:))) &
        + 0.5_dp * (1.0_dp - alpha) * sum(x(first_penalized:)**2)
    end if
  end function regularizer_en

end module rpeglmen_penalty
