! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_wrappers
  use nlme_kinds, only : dp
  use segmented_status
  use segmented_types
  use segmented_utils, only : quantile_value
  use segmented_fit
  use segmented_mixed, only : fit_segmented_lme
  implicit none
  private
  public :: segmented_lm, segmented_glm, stepmented_lm, stepmented_glm
  public :: segmented_lme, segmented_numeric, stepmented_numeric, segreg, stepreg
contains
  subroutine segmented_lm(y, x, z, psi0, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_segmented_lm(y, x, z, psi0, result, weights, offset, control)
  end subroutine segmented_lm

  subroutine segmented_glm(y, x, z, psi0, family, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_segmented_glm(y, x, z, psi0, family, result, weights, offset, control)
  end subroutine segmented_glm

  subroutine stepmented_lm(y, x, z, psi0, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_stepmented_lm(y, x, z, psi0, result, weights, offset, control)
  end subroutine stepmented_lm

  subroutine stepmented_glm(y, x, z, psi0, family, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_stepmented_glm(y, x, z, psi0, family, result, weights, offset, control)
  end subroutine stepmented_glm

  subroutine segmented_lme(y, x, z, psi0, random_design, group, result, options, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:), random_design(:,:)
    integer, intent(in) :: group(:)
    type(segmented_lme_result), intent(out) :: result
    type(segmented_lme_options), intent(in), optional :: options
    type(segmented_control), intent(in), optional :: control
    call fit_segmented_lme(y, x, z, psi0, random_design, group, result, options, control)
  end subroutine segmented_lme

  subroutine segmented_numeric(y, result, x_values, psi0, n_break, weights, control)
    real(dp), intent(in) :: y(:)
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: x_values(:), psi0(:), weights(:)
    integer, intent(in), optional :: n_break
    type(segmented_control), intent(in), optional :: control
    call numeric_fit(y, SEGMENTED_CONTINUOUS, result, x_values, psi0, n_break, &
        weights, control)
  end subroutine segmented_numeric

  subroutine stepmented_numeric(y, result, x_values, psi0, n_break, weights, control)
    real(dp), intent(in) :: y(:)
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: x_values(:), psi0(:), weights(:)
    integer, intent(in), optional :: n_break
    type(segmented_control), intent(in), optional :: control
    call numeric_fit(y, SEGMENTED_STEP, result, x_values, psi0, n_break, weights, control)
  end subroutine stepmented_numeric

  subroutine segreg(y, x, z, psi0, family, kind, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family, kind
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    if (kind == SEGMENTED_CONTINUOUS) then
      call fit_segmented_glm(y, x, z, psi0, family, result, weights, offset, control)
    else if (kind == SEGMENTED_STEP) then
      call fit_stepmented_glm(y, x, z, psi0, family, result, weights, offset, control)
    else
      result%status = SEG_INVALID_ARGUMENT
    end if
  end subroutine segreg

  subroutine stepreg(y, x, z, psi0, family, result, weights, offset, control)
    real(dp), intent(in) :: y(:), x(:,:), z(:,:), psi0(:)
    integer, intent(in) :: family
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:), offset(:)
    type(segmented_control), intent(in), optional :: control
    call fit_stepmented_glm(y, x, z, psi0, family, result, weights, offset, control)
  end subroutine stepreg

  subroutine numeric_fit(y, kind, result, x_values, psi0, n_break, weights, control)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: kind
    type(segmented_result), intent(out) :: result
    real(dp), intent(in), optional :: x_values(:), psi0(:), weights(:)
    integer, intent(in), optional :: n_break
    type(segmented_control), intent(in), optional :: control
    real(dp), allocatable :: xval(:), x(:,:), z(:,:), starts(:)
    integer :: n, m, i, j
    n = size(y)
    if (n < 5) then
      result%status = SEG_DIMENSION_ERROR
      return
    end if
    allocate(xval(n))
    if (present(x_values)) then
      if (size(x_values) /= n) then
        result%status = SEG_DIMENSION_ERROR
        return
      end if
      xval = x_values
    else
      do i = 1, n
        xval(i) = real(i, dp) / real(n, dp)
      end do
    end if
    if (present(psi0)) then
      m = size(psi0)
    else
      m = 1
      if (present(n_break)) m = n_break
    end if
    if (m < 1) then
      result%status = SEG_INVALID_ARGUMENT
      return
    end if
    allocate(x(n, 2), z(n, m), starts(m))
    x(:, 1) = 1.0_dp
    x(:, 2) = xval
    do j = 1, m
      z(:, j) = xval
      starts(j) = quantile_value(xval, real(j, dp) / real(m + 1, dp))
    end do
    if (present(psi0)) starts = psi0
    if (kind == SEGMENTED_CONTINUOUS) then
      call fit_segmented_lm(y, x, z, starts, result, weights=weights, control=control)
    else
      call fit_stepmented_lm(y, x(:, :1), z, starts, result, weights=weights, control=control)
    end if
  end subroutine numeric_fit
end module segmented_wrappers
