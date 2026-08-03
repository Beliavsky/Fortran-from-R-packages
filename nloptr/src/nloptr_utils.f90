! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_utils
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nloptr_kinds, only: dp
  use nloptr_types, only: nloptr_problem, nloptr_options, NLOPT_INVALID_ARGS
  implicit none
  private
  public :: project_bounds, validate_problem, max_abs, vec_norm, outer_product
  public :: halton_point, finite_vector, constraint_violation, status_message

contains

  subroutine project_bounds(x, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    integer :: i
    do i = 1, size(x)
      if (ieee_is_finite(lower(i))) x(i) = max(x(i), lower(i))
      if (ieee_is_finite(upper(i))) x(i) = min(x(i), upper(i))
    end do
  end subroutine project_bounds

  integer function validate_problem(problem, x0) result(status)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    status = 0
    if (.not. associated(problem%objective)) then
      status = NLOPT_INVALID_ARGS
      return
    end if
    if (problem%n <= 0 .or. size(x0) /= problem%n) then
      status = NLOPT_INVALID_ARGS
      return
    end if
    if (.not. allocated(problem%lower) .or. .not. allocated(problem%upper)) then
      status = NLOPT_INVALID_ARGS
      return
    end if
    if (size(problem%lower) /= problem%n .or. size(problem%upper) /= problem%n) then
      status = NLOPT_INVALID_ARGS
      return
    end if
    if (any(problem%lower > problem%upper) .or. .not. finite_vector(x0)) then
      status = NLOPT_INVALID_ARGS
      return
    end if
    if (problem%n_ineq > 0 .and. .not. associated(problem%inequality)) status = NLOPT_INVALID_ARGS
    if (problem%n_eq > 0 .and. .not. associated(problem%equality)) status = NLOPT_INVALID_ARGS
  end function validate_problem

  pure real(dp) function max_abs(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = maxval(abs(x))
    end if
  end function max_abs

  pure real(dp) function vec_norm(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vec_norm

  pure function outer_product(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))
    integer :: j
    do j = 1, size(b)
      c(:, j) = a * b(j)
    end do
  end function outer_product

  pure logical function finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    ok = .true.
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) then
        ok = .false.
        return
      end if
    end do
  end function finite_vector

  pure real(dp) function radical_inverse(index, base) result(value)
    integer, intent(in) :: index, base
    integer :: i
    real(dp) :: factor
    i = max(0, index)
    value = 0.0_dp
    factor = 1.0_dp / real(base, dp)
    do while (i > 0)
      value = value + factor * real(mod(i, base), dp)
      i = i / base
      factor = factor / real(base, dp)
    end do
  end function radical_inverse

  subroutine halton_point(index, point)
    integer, intent(in) :: index
    real(dp), intent(out) :: point(:)
    integer, parameter :: primes(20) = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, &
      31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    integer :: j
    do j = 1, size(point)
      if (j <= size(primes)) then
        point(j) = radical_inverse(index, primes(j))
      else
        point(j) = radical_inverse(index, 73 + 2 * (j - size(primes) - 1))
      end if
    end do
  end subroutine halton_point

  pure real(dp) function constraint_violation(g, h) result(value)
    real(dp), intent(in) :: g(:), h(:)
    value = 0.0_dp
    if (size(g) > 0) value = max(value, maxval(max(0.0_dp, g)))
    if (size(h) > 0) value = max(value, maxval(abs(h)))
  end function constraint_violation

  pure function status_message(status) result(message)
    integer, intent(in) :: status
    character(len=160) :: message
    select case (status)
    case (1)
      message = 'success'
    case (2)
      message = 'stop value reached'
    case (3)
      message = 'function tolerance reached'
    case (4)
      message = 'parameter tolerance reached'
    case (5)
      message = 'maximum evaluations reached'
    case (6)
      message = 'maximum time reached'
    case (-2)
      message = 'invalid arguments'
    case (-4)
      message = 'roundoff limited progress'
    case (-5)
      message = 'forced stop requested by callback'
    case default
      message = 'optimization failure'
    end select
  end function status_message
end module nloptr_utils
