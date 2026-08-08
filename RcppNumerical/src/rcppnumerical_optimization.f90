module rcppnumerical_optimization
  use rcppnumerical_kinds, only : dp
  use rcppnumerical_callbacks, only : objective_gradient_interface
  use lbfgs, only : lbfgs_parameter_t, lbfgs_result_t, lbfgs_minimize, &
                    lbfgs_linesearch_backtracking_strong_wolfe
  use lbfgsb3_mod, only : lbfgsb_control_t, lbfgsb_result_t, lbfgsb_minimize
  implicit none
  private

  type, public :: optimization_result_t
    real(dp) :: value = huge(1.0_dp)
    integer :: status = -1
    integer :: native_status = -1
    integer :: iterations = 0
    integer :: evaluations = 0
    real(dp) :: gradient_norm = huge(1.0_dp)
    logical :: converged = .false.
  end type optimization_result_t

  public :: optim_lbfgs, optim_lbfgsb

contains

  subroutine optim_lbfgs(f, x, result, maxit, eps_f, eps_g, user_data)
    procedure(objective_gradient_interface) :: f
    real(dp), intent(inout) :: x(:)
    type(optimization_result_t), intent(out) :: result
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: eps_f, eps_g
    class(*), intent(inout), optional :: user_data

    type(lbfgs_parameter_t) :: param
    type(lbfgs_result_t) :: native

    param = lbfgs_parameter_t()
    param%epsilon = 1.0e-5_dp
    if (present(eps_g)) param%epsilon = eps_g
    param%past = 1
    param%delta = 1.0e-6_dp
    if (present(eps_f)) param%delta = eps_f
    param%max_iterations = 300
    if (present(maxit)) param%max_iterations = maxit
    param%max_linesearch = 100
    param%linesearch = lbfgs_linesearch_backtracking_strong_wolfe

    if (present(user_data)) then
      call lbfgs_minimize(adapter, x, native, param, user_data=user_data)
    else
      call lbfgs_minimize(adapter, x, native, param)
    end if
    result%value = native%value
    result%native_status = native%status
    result%iterations = native%iterations
    result%evaluations = native%evaluations
    result%gradient_norm = native%gradient_norm
    result%converged = native%status >= 0
    if (result%converged) then
      result%status = 0
    else
      result%status = -1
    end if

  contains

    subroutine adapter(x_current, value, gradient, step, callback_data)
      real(dp), intent(in) :: x_current(:)
      real(dp), intent(out) :: value
      real(dp), intent(out) :: gradient(:)
      real(dp), intent(in) :: step
      class(*), intent(inout), optional :: callback_data
      if (present(callback_data)) then
        call f(x_current, value, gradient, callback_data)
      else
        call f(x_current, value, gradient)
      end if
    end subroutine adapter

  end subroutine optim_lbfgs

  subroutine optim_lbfgsb(f, x, lower, upper, result, maxit, eps_f, eps_g, user_data)
    procedure(objective_gradient_interface) :: f
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:), upper(:)
    type(optimization_result_t), intent(out) :: result
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: eps_f, eps_g
    class(*), intent(inout), optional :: user_data

    type(lbfgsb_control_t) :: control
    type(lbfgsb_result_t) :: native
    real(dp) :: rel_f

    control = lbfgsb_control_t()
    control%memory = 6
    control%pgtol = 1.0e-5_dp
    if (present(eps_g)) control%pgtol = eps_g
    rel_f = 1.0e-6_dp
    if (present(eps_f)) rel_f = eps_f
    control%factr = max(0.0_dp, rel_f/epsilon(1.0_dp))
    control%max_iterations = 300
    if (present(maxit)) control%max_iterations = maxit
    control%max_evaluations = max(1000, 100*(control%max_iterations + 1))
    control%abstol = 0.0_dp
    control%reltol = 0.0_dp

    if (present(user_data)) then
      call lbfgsb_minimize(x, adapter, native, lower, upper, control, user_data)
    else
      call lbfgsb_minimize(x, adapter, native, lower, upper, control)
    end if
    result%value = native%value
    result%native_status = native%task
    result%iterations = native%iterations
    result%evaluations = native%function_evaluations
    result%gradient_norm = native%projected_gradient_norm
    result%converged = native%success
    if (result%converged) then
      result%status = 0
    else
      result%status = -1
    end if

  contains

    subroutine adapter(x_current, value, gradient, callback_data)
      real(dp), intent(in) :: x_current(:)
      real(dp), intent(out) :: value
      real(dp), intent(out) :: gradient(:)
      class(*), intent(inout), optional :: callback_data
      if (present(callback_data)) then
        call f(x_current, value, gradient, callback_data)
      else
        call f(x_current, value, gradient)
      end if
    end subroutine adapter

  end subroutine optim_lbfgsb

end module rcppnumerical_optimization
