! SPDX-License-Identifier: GPL-3.0-only
module rgenoud
  use rgenoud_kinds, only : dp
  use rgenoud_types, only : genoud_options, genoud_result, objective_fn, gradient_fn, &
    lexical_objective_fn, lexical_better_fn
  use rgenoud_derivatives, only : numerical_gradient, numerical_hessian
  use rgenoud_core, only : run_genoud_core
  use rgenoud_stats, only : sample_moments
  implicit none
  private
  public :: dp, genoud_options, genoud_result
  public :: genoud_optimize, genoud_optimize_lexical
  public :: numerical_gradient, numerical_hessian, sample_moments
contains
  subroutine dispatch_gradient(callback, x, g)
    procedure(gradient_fn) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)

    call callback(x, g)
  end subroutine dispatch_gradient

  real(dp) function dispatch_objective(callback, x) result(value)
    procedure(objective_fn) :: callback
    real(dp), intent(in) :: x(:)

    value = callback(x)
  end function dispatch_objective

  logical function dispatch_comparator(callback, a, b) result(is_better)
    procedure(lexical_better_fn) :: callback
    real(dp), intent(in) :: a(:), b(:)

    is_better = callback(a, b)
  end function dispatch_comparator

  subroutine dispatch_lexical_objective(callback, x, f)
    procedure(lexical_objective_fn) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)

    call callback(x, f)
  end subroutine dispatch_lexical_objective

  subroutine genoud_optimize(fn, lower, upper, options, result, gradient, &
      starting_values, hessian)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: lower(:), upper(:)
    type(genoud_options), intent(in), optional :: options
    type(genoud_result), intent(out) :: result
    procedure(gradient_fn), optional :: gradient
    real(dp), intent(in), optional :: starting_values(:, :)
    logical, intent(in), optional :: hessian
    type(genoud_options) :: opt
    logical :: want_hessian

    opt = genoud_options()
    if (present(options)) opt = options
    want_hessian = .false.
    if (present(hessian)) want_hessian = hessian

    if (present(starting_values)) then
      call run_genoud_core(eval_adapter, fn, grad_adapter, better_adapter, lower, upper, &
        1, opt, result, .not. opt%integer_parameters, .not. opt%integer_parameters, &
        starting_values)
    else
      call run_genoud_core(eval_adapter, fn, grad_adapter, better_adapter, lower, upper, &
        1, opt, result, .not. opt%integer_parameters, .not. opt%integer_parameters)
    end if
    if (want_hessian .and. .not. opt%integer_parameters) then
      allocate(result%hessian(size(lower), size(lower)))
      call numerical_hessian(fn, result%par, result%hessian, lower, upper, &
        opt%boundary_enforcement == 2)
    end if

  contains
    subroutine eval_adapter(x, f)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
      f(1) = dispatch_objective(fn, x)
      if (.not. ieee_finite_scalar(f(1))) then
        f(1) = merge(-huge(1.0_dp), huge(1.0_dp), opt%maximize)
      end if
    end subroutine eval_adapter

    subroutine grad_adapter(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      if (present(gradient)) then
        call dispatch_gradient(gradient, x, g)
      else
        call numerical_gradient(fn, x, g, lower, upper, opt%boundary_enforcement == 2)
      end if
    end subroutine grad_adapter

    logical function better_adapter(a, b) result(is_better)
      real(dp), intent(in) :: a(:), b(:)
      if (opt%maximize) then
        is_better = a(1) > b(1)
      else
        is_better = a(1) < b(1)
      end if
    end function better_adapter

    logical function ieee_finite_scalar(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
    end function ieee_finite_scalar
  end subroutine genoud_optimize

  subroutine genoud_optimize_lexical(fn, nfit, lower, upper, options, result, &
      comparator, local_objective, local_gradient, starting_values, hessian)
    procedure(lexical_objective_fn) :: fn
    integer, intent(in) :: nfit
    real(dp), intent(in) :: lower(:), upper(:)
    type(genoud_options), intent(in), optional :: options
    type(genoud_result), intent(out) :: result
    procedure(lexical_better_fn), optional :: comparator
    procedure(objective_fn), optional :: local_objective
    procedure(gradient_fn), optional :: local_gradient
    real(dp), intent(in), optional :: starting_values(:, :)
    logical, intent(in), optional :: hessian
    type(genoud_options) :: opt
    logical :: local_enabled, want_hessian

    opt = genoud_options()
    if (present(options)) opt = options
    local_enabled = present(local_objective) .and. .not. opt%integer_parameters
    if (.not. local_enabled) then
      opt%use_bfgs = .false.
      opt%gradient_check = .false.
      opt%operator_weights(9) = 0.0_dp
    end if
    want_hessian = .false.
    if (present(hessian)) want_hessian = hessian

    if (present(starting_values)) then
      call run_genoud_core(fn, local_adapter, grad_adapter, better_adapter, lower, upper, &
        nfit, opt, result, local_enabled, local_enabled, starting_values)
    else
      call run_genoud_core(fn, local_adapter, grad_adapter, better_adapter, lower, upper, &
        nfit, opt, result, local_enabled, local_enabled)
    end if

    if (want_hessian .and. local_enabled) then
      allocate(result%hessian(size(lower), size(lower)))
      call numerical_hessian(local_adapter, result%par, result%hessian, lower, upper, &
        opt%boundary_enforcement == 2)
    end if

  contains
    real(dp) function local_adapter(x) result(fval)
      real(dp), intent(in) :: x(:)
      real(dp) :: fv(nfit)
      if (present(local_objective)) then
        fval = dispatch_objective(local_objective, x)
      else
        call dispatch_lexical_objective(fn, x, fv)
        fval = fv(1)
      end if
    end function local_adapter

    subroutine grad_adapter(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      if (present(local_gradient)) then
        call dispatch_gradient(local_gradient, x, g)
      else
        call numerical_gradient(local_adapter, x, g, lower, upper, &
          opt%boundary_enforcement == 2)
      end if
    end subroutine grad_adapter

    logical function better_adapter(a, b) result(is_better)
      real(dp), intent(in) :: a(:), b(:)
      integer :: i
      if (present(comparator)) then
        is_better = dispatch_comparator(comparator, a, b)
        return
      end if
      is_better = .false.
      do i = 1, min(size(a), size(b))
        if (abs(a(i) - b(i)) <= 0.0_dp) cycle
        if (opt%maximize) then
          is_better = a(i) > b(i)
        else
          is_better = a(i) < b(i)
        end if
        return
      end do
    end function better_adapter
  end subroutine genoud_optimize_lexical
end module rgenoud
