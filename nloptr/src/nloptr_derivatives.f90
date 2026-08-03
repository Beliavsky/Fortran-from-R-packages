! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_derivatives
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nloptr_kinds, only: dp
  use nloptr_types, only: objective_callback, vector_callback, derivative_check_result
  implicit none
  private
  public :: nl_grad, nl_jacobian, check_derivatives

contains

  subroutine nl_grad(x, func, gradient, status, step)
    real(dp), intent(in) :: x(:)
    procedure(objective_callback) :: func
    real(dp), intent(out) :: gradient(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    real(dp) :: xp(size(x)), xm(size(x)), dummy(size(x)), fp, fm, h
    integer :: i, statp, statm

    status = 0
    if (size(gradient) /= size(x) .or. size(x) == 0) then
      status = -2
      return
    end if
    do i = 1, size(x)
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(i)))
      if (present(step)) h = max(abs(step), epsilon(1.0_dp)) * max(1.0_dp, abs(x(i)))
      xp = x
      xm = x
      xp(i) = xp(i) + h
      xm(i) = xm(i) - h
      dummy = 0.0_dp
      call func(xp, fp, dummy, .false., statp)
      dummy = 0.0_dp
      call func(xm, fm, dummy, .false., statm)
      if (statp /= 0 .or. statm /= 0 .or. .not. ieee_is_finite(fp) .or. .not. ieee_is_finite(fm)) then
        status = -5
        return
      end if
      gradient(i) = (fp - fm) / (2.0_dp * h)
    end do
  end subroutine nl_grad

  subroutine nl_jacobian(x, m, func, jacobian, status, step)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: m
    procedure(vector_callback) :: func
    real(dp), intent(out) :: jacobian(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    real(dp) :: xp(size(x)), xm(size(x)), fp(m), fm(m), dummy(m, size(x)), h
    integer :: i, statp, statm

    status = 0
    if (m < 0 .or. size(jacobian, 1) /= m .or. size(jacobian, 2) /= size(x)) then
      status = -2
      return
    end if
    do i = 1, size(x)
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(i)))
      if (present(step)) h = max(abs(step), epsilon(1.0_dp)) * max(1.0_dp, abs(x(i)))
      xp = x
      xm = x
      xp(i) = xp(i) + h
      xm(i) = xm(i) - h
      dummy = 0.0_dp
      call func(xp, fp, dummy, .false., statp)
      dummy = 0.0_dp
      call func(xm, fm, dummy, .false., statm)
      if (statp /= 0 .or. statm /= 0 .or. any(.not. ieee_is_finite(fp)) .or. &
          any(.not. ieee_is_finite(fm))) then
        status = -5
        return
      end if
      jacobian(:, i) = (fp - fm) / (2.0_dp * h)
    end do
  end subroutine nl_jacobian

  subroutine check_derivatives(x, func, analytic, result, tolerance)
    real(dp), intent(in) :: x(:)
    procedure(objective_callback) :: func
    real(dp), intent(in) :: analytic(:)
    type(derivative_check_result), intent(out) :: result
    real(dp), intent(in), optional :: tolerance
    real(dp) :: tol, denom
    integer :: i, status

    tol = 1.0e-4_dp
    if (present(tolerance)) tol = max(tolerance, 0.0_dp)
    allocate(result%analytic(1, size(x)), result%numeric(1, size(x)))
    allocate(result%relative_error(1, size(x)), result%warning(1, size(x)))
    result%analytic(1, :) = analytic
    call nl_grad(x, func, result%numeric(1, :), status)
    result%status = status
    result%n_warnings = 0
    if (status /= 0 .or. size(analytic) /= size(x)) then
      result%warning = .true.
      result%n_warnings = size(x)
      return
    end if
    do i = 1, size(x)
      denom = max(1.0_dp, abs(analytic(i)), abs(result%numeric(1, i)))
      result%relative_error(1, i) = abs(analytic(i) - result%numeric(1, i)) / denom
      result%warning(1, i) = result%relative_error(1, i) > tol
      if (result%warning(1, i)) result%n_warnings = result%n_warnings + 1
    end do
  end subroutine check_derivatives
end module nloptr_derivatives
