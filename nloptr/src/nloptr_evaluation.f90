! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_evaluation
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use nloptr_kinds, only: dp
  use nloptr_types, only: nloptr_problem
  use nloptr_derivatives, only: nl_grad, nl_jacobian
  use nloptr_utils, only: constraint_violation
  implicit none
  private
  public :: evaluate_problem, evaluate_penalty

contains

  subroutine evaluate_problem(problem, x, f, grad, g, jg, h, jh, need_gradient, status, evals)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: grad(:)
    real(dp), intent(out) :: g(:), jg(:, :), h(:), jh(:, :)
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    integer, intent(inout) :: evals
    real(dp) :: test_grad(size(x))
    integer :: stat

    status = 0
    if (need_gradient) then
      grad = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(jg) > 0) jg = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(jh) > 0) jh = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      grad = 0.0_dp
      if (size(jg) > 0) jg = 0.0_dp
      if (size(jh) > 0) jh = 0.0_dp
    end if
    if (size(g) > 0) g = 0.0_dp
    if (size(h) > 0) h = 0.0_dp

    call problem%objective(x, f, grad, need_gradient, stat)
    evals = evals + 1
    if (stat /= 0 .or. .not. ieee_is_finite(f)) then
      status = merge(stat, -5, stat /= 0)
      return
    end if
    if (need_gradient) then
      if (any(.not. ieee_is_finite(grad))) then
        call nl_grad(x, problem%objective, test_grad, stat)
        if (stat /= 0) then
          status = stat
          return
        end if
        grad = test_grad
        evals = evals + 2 * size(x)
      end if
    end if

    if (problem%n_ineq > 0) then
      call problem%inequality(x, g, jg, need_gradient, stat)
      evals = evals + 1
      if (stat /= 0 .or. any(.not. ieee_is_finite(g))) then
        status = merge(stat, -5, stat /= 0)
        return
      end if
      if (need_gradient .and. any(.not. ieee_is_finite(jg))) then
        call nl_jacobian(x, problem%n_ineq, problem%inequality, jg, stat)
        if (stat /= 0) then
          status = stat
          return
        end if
        evals = evals + 2 * size(x)
      end if
    end if

    if (problem%n_eq > 0) then
      call problem%equality(x, h, jh, need_gradient, stat)
      evals = evals + 1
      if (stat /= 0 .or. any(.not. ieee_is_finite(h))) then
        status = merge(stat, -5, stat /= 0)
        return
      end if
      if (need_gradient .and. any(.not. ieee_is_finite(jh))) then
        call nl_jacobian(x, problem%n_eq, problem%equality, jh, stat)
        if (stat /= 0) then
          status = stat
          return
        end if
        evals = evals + 2 * size(x)
      end if
    end if
  end subroutine evaluate_problem

  subroutine evaluate_penalty(problem, x, rho, value, gradient, violation, need_gradient, status, evals)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:), rho
    real(dp), intent(out) :: value
    real(dp), intent(inout) :: gradient(:), violation
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    integer, intent(inout) :: evals
    real(dp) :: f, grad(size(x))
    real(dp), allocatable :: g(:), h(:), jg(:, :), jh(:, :), gp(:)

    allocate(g(problem%n_ineq), h(problem%n_eq))
    allocate(jg(problem%n_ineq, problem%n), jh(problem%n_eq, problem%n))
    call evaluate_problem(problem, x, f, grad, g, jg, h, jh, need_gradient, status, evals)
    if (status /= 0) return
    violation = constraint_violation(g, h)
    value = f
    gradient = grad
    if (problem%n_ineq > 0) then
      allocate(gp(problem%n_ineq))
      gp = max(0.0_dp, g)
      value = value + rho * dot_product(gp, gp)
      if (need_gradient) gradient = gradient + 2.0_dp * rho * matmul(transpose(jg), gp)
    end if
    if (problem%n_eq > 0) then
      value = value + rho * dot_product(h, h)
      if (need_gradient) gradient = gradient + 2.0_dp * rho * matmul(transpose(jh), h)
    end if
  end subroutine evaluate_penalty
end module nloptr_evaluation
