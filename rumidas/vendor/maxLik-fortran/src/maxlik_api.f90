! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_api
  use maxlik_kinds, only: dp
  use maxlik_types, only: maxlik_problem, maxlik_control, maxlik_result
  use maxlik_status
  use maxlik_solvers
  use maxlik_evaluation, only: constraint_violation
  use maxlik_inference, only: finalize_result
  implicit none
  private

  public :: max_lik

contains

  subroutine max_lik(problem, start, result, method, control)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    type(maxlik_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    type(maxlik_control), intent(in), optional :: control

    type(maxlik_control) :: ctrl
    type(maxlik_problem) :: work, raw_problem
    type(maxlik_result) :: trial = maxlik_result()
    real(dp), allocatable :: x(:)
    real(dp) :: rho, violation
    integer :: outer, total_function, total_gradient, total_hessian
    character(len=32) :: selected
    logical :: constrained

    ctrl = maxlik_control()
    if (present(control)) ctrl = control
    selected = 'nr'
    if (present(method)) selected = lower_ascii(trim(method))
    constrained = allocated(problem%eq_a) .or. allocated(problem%ineq_a)

    if (.not. constrained) then
      call dispatch_solver(problem, start, selected, ctrl, result)
      raw_problem = problem
      raw_problem%penalty_rho = 0.0_dp
      call finalize_result(raw_problem, ctrl, result)
      return
    end if

    allocate(x(size(start)))
    x = start
    rho = max(ctrl%constraint_rho0, tiny(1.0_dp))
    total_function = 0
    total_gradient = 0
    total_hessian = 0
    violation = huge(1.0_dp)
    do outer = 1, ctrl%constraint_max_outer
      work = problem
      work%penalty_rho = rho
      call dispatch_solver(work, x, selected, ctrl, trial)
      total_function = total_function + trial%function_count
      total_gradient = total_gradient + trial%gradient_count
      total_hessian = total_hessian + trial%hessian_count
      if (.not. allocated(trial%estimate)) then
        result = trial
        return
      end if
      x = trial%estimate
      violation = constraint_violation(problem, x)
      result = trial
      result%outer_iterations = outer
      result%constraint_violation = violation
      if (violation <= ctrl%constraint_tol) exit
      rho = rho * ctrl%constraint_rho_factor
    end do
    result%function_count = total_function
    result%gradient_count = total_gradient
    result%hessian_count = total_hessian
    if (violation > ctrl%constraint_tol) then
      result%code = MAXLIK_CONSTRAINT_FAILURE
      result%message = maxlik_message(result%code)
      result%converged = .false.
    end if
    raw_problem = problem
    raw_problem%penalty_rho = 0.0_dp
    call finalize_result(raw_problem, ctrl, result)
    result%constraint_violation = violation
  end subroutine max_lik

  subroutine dispatch_solver(problem, start, method, control, result)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: start(:)
    character(len=*), intent(in) :: method
    type(maxlik_control), intent(in) :: control
    type(maxlik_result), intent(out) :: result

    select case (trim(lower_ascii(method)))
    case ('nr', 'newton', 'newton-raphson')
      call solve_newton(problem, start, control, result)
    case ('bfgs')
      call solve_bfgs(problem, start, control, result)
    case ('bfgsr', 'bfgs-r')
      call solve_bfgsr(problem, start, control, result)
    case ('bhhh')
      call solve_bhhh(problem, start, control, result)
    case ('cg', 'conjugate-gradient')
      call solve_cg(problem, start, control, result)
    case ('nm', 'nelder-mead', 'nelder_mead')
      call solve_nelder_mead(problem, start, control, result)
    case ('sann', 'simulated-annealing')
      call solve_sann(problem, start, control, result)
    case ('sga', 'stochastic-gradient')
      call solve_sga(problem, start, control, result)
    case ('adam')
      call solve_adam(problem, start, control, result)
    case default
      result%code = MAXLIK_INVALID_INPUT
      result%message = 'unknown maximization method: ' // trim(method)
      result%method = trim(method)
    end select
  end subroutine dispatch_solver

  pure function lower_ascii(text) result(lower)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lower
    integer :: i, code

    lower = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
    end do
  end function lower_ascii

end module maxlik_api
