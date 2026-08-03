! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_api
  use nloptr_kinds, only: dp
  use nloptr_types
  use nloptr_solvers, only: optimize_problem
  use nloptr_derivatives, only: nl_grad, nl_jacobian, check_derivatives
  implicit none
  private

  public :: dp, nloptr_problem, nloptr_options, nloptr_result, derivative_check_result
  public :: objective_callback, vector_callback
  public :: nloptr, lbfgs, varmetric, tnewton, neldermead, sbplx
  public :: cobyla, bobyqa, newuoa, slsqp, mma, ccsaq, auglag
  public :: direct, direct_l, crs2lm, isres, stogo, mlsl
  public :: nl_grad, nl_jacobian, check_derivatives
  public :: nl_opts, nloptr_get_default_options, nloptr_print_options, is_nloptr
  public :: NLOPT_FAILURE, NLOPT_INVALID_ARGS, NLOPT_OUT_OF_MEMORY
  public :: NLOPT_ROUNDOFF_LIMITED, NLOPT_FORCED_STOP, NLOPT_SUCCESS
  public :: NLOPT_STOPVAL_REACHED, NLOPT_FTOL_REACHED, NLOPT_XTOL_REACHED
  public :: NLOPT_MAXEVAL_REACHED, NLOPT_MAXTIME_REACHED

contains

  subroutine nloptr(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    type(nloptr_options) :: actual
    actual = nl_opts()
    if (present(options)) actual = options
    call optimize_problem(problem, x0, actual, result)
  end subroutine nloptr

  function nl_opts(algorithm) result(options)
    character(len=*), intent(in), optional :: algorithm
    type(nloptr_options) :: options
    if (present(algorithm)) options%algorithm = algorithm
  end function nl_opts

  function nloptr_get_default_options() result(options)
    type(nloptr_options) :: options
    options = nloptr_options()
  end function nloptr_get_default_options

  subroutine nloptr_print_options(options, unit)
    type(nloptr_options), intent(in), optional :: options
    integer, intent(in), optional :: unit
    type(nloptr_options) :: actual
    integer :: output_unit
    actual = nl_opts()
    if (present(options)) actual = options
    output_unit = 6
    if (present(unit)) output_unit = unit
    write(output_unit, '(a,a)') 'algorithm: ', trim(actual%algorithm)
    write(output_unit, '(a,es12.4)') 'ftol_rel: ', actual%ftol_rel
    write(output_unit, '(a,es12.4)') 'xtol_rel: ', actual%xtol_rel
    write(output_unit, '(a,i0)') 'maxeval: ', actual%maxeval
    write(output_unit, '(a,es12.4)') 'constraint_tol: ', actual%constraint_tol
    write(output_unit, '(a,i0)') 'population: ', actual%population
    write(output_unit, '(a,i0)') 'seed: ', actual%seed
  end subroutine nloptr_print_options

  pure logical function is_nloptr(result) result(answer)
    type(nloptr_result), intent(in) :: result
    answer = allocated(result%solution) .and. result%status /= NLOPT_FAILURE
  end function is_nloptr

  subroutine run_named(problem, x0, name, options, result, local_name)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    character(len=*), intent(in) :: name
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    character(len=*), intent(in), optional :: local_name
    type(nloptr_options) :: actual
    actual = nl_opts(name)
    if (present(options)) then
      actual = options
      actual%algorithm = name
    end if
    if (present(local_name)) actual%local_algorithm = local_name
    call optimize_problem(problem, x0, actual, result)
  end subroutine run_named

  subroutine lbfgs(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LD_LBFGS', options, result)
  end subroutine lbfgs

  subroutine varmetric(problem, x0, options, result, rank_one)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    logical, intent(in), optional :: rank_one
    character(len=40) :: name
    name = 'NLOPT_LD_VAR2'
    if (present(rank_one)) then
      if (rank_one) name = 'NLOPT_LD_VAR1'
    end if
    call run_named(problem, x0, name, options, result)
  end subroutine varmetric

  subroutine tnewton(problem, x0, options, result, restart, precondition)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    logical, intent(in), optional :: restart, precondition
    logical :: use_restart, use_precondition
    character(len=40) :: name
    use_restart = .false.
    use_precondition = .false.
    if (present(restart)) use_restart = restart
    if (present(precondition)) use_precondition = precondition
    if (use_restart .and. use_precondition) then
      name = 'NLOPT_LD_TNEWTON_PRECOND_RESTART'
    else if (use_restart) then
      name = 'NLOPT_LD_TNEWTON_RESTART'
    else if (use_precondition) then
      name = 'NLOPT_LD_TNEWTON_PRECOND'
    else
      name = 'NLOPT_LD_TNEWTON'
    end if
    call run_named(problem, x0, name, options, result)
  end subroutine tnewton

  subroutine neldermead(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LN_NELDERMEAD', options, result)
  end subroutine neldermead

  subroutine sbplx(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LN_SBPLX', options, result)
  end subroutine sbplx

  subroutine cobyla(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LN_COBYLA', options, result)
  end subroutine cobyla

  subroutine bobyqa(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LN_BOBYQA', options, result)
  end subroutine bobyqa

  subroutine newuoa(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LN_NEWUOA', options, result)
  end subroutine newuoa

  subroutine slsqp(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LD_SLSQP', options, result)
  end subroutine slsqp

  subroutine mma(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LD_MMA', options, result)
  end subroutine mma

  subroutine ccsaq(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_LD_CCSAQ', options, result)
  end subroutine ccsaq

  subroutine auglag(problem, x0, options, result, derivative_free)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    logical, intent(in), optional :: derivative_free
    character(len=40) :: name, local_name
    name = 'NLOPT_LD_AUGLAG'
    local_name = 'NLOPT_LD_LBFGS'
    if (present(derivative_free)) then
      if (derivative_free) then
        name = 'NLOPT_LN_AUGLAG'
        local_name = 'NLOPT_LN_NELDERMEAD'
      end if
    end if
    call run_named(problem, x0, name, options, result, local_name)
  end subroutine auglag

  subroutine direct(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_GN_DIRECT', options, result, 'NLOPT_LN_NELDERMEAD')
  end subroutine direct

  subroutine direct_l(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_GN_DIRECT_L', options, result, 'NLOPT_LN_NELDERMEAD')
  end subroutine direct_l

  subroutine crs2lm(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_GN_CRS2_LM', options, result, 'NLOPT_LN_NELDERMEAD')
  end subroutine crs2lm

  subroutine isres(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    call run_named(problem, x0, 'NLOPT_GN_ISRES', options, result, 'NLOPT_LN_NELDERMEAD')
  end subroutine isres

  subroutine stogo(problem, x0, options, result, randomized)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    logical, intent(in), optional :: randomized
    character(len=40) :: name
    name = 'NLOPT_GD_STOGO'
    if (present(randomized)) then
      if (randomized) name = 'NLOPT_GD_STOGO_RAND'
    end if
    call run_named(problem, x0, name, options, result, 'NLOPT_LD_LBFGS')
  end subroutine stogo

  subroutine mlsl(problem, x0, options, result, derivative_free)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in), optional :: options
    type(nloptr_result), intent(out) :: result
    logical, intent(in), optional :: derivative_free
    character(len=40) :: name, local_name
    name = 'NLOPT_GD_MLSL_LDS'
    local_name = 'NLOPT_LD_LBFGS'
    if (present(derivative_free)) then
      if (derivative_free) then
        name = 'NLOPT_GN_MLSL_LDS'
        local_name = 'NLOPT_LN_NELDERMEAD'
      end if
    end if
    call run_named(problem, x0, name, options, result, local_name)
  end subroutine mlsl
end module nloptr_api
