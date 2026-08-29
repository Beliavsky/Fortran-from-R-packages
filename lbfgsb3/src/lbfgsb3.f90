module lbfgsb3_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
       ieee_positive_inf, ieee_negative_inf
  use, intrinsic :: iso_fortran_env, only : output_unit
  use lbfgsb3_core_mod, only : setulb, set_core_output
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  integer, parameter, public :: lbfgsb_success = 0
  integer, parameter, public :: lbfgsb_iteration_limit = 1
  integer, parameter, public :: lbfgsb_warning = 51
  integer, parameter, public :: lbfgsb_input_error = 52
  integer, parameter, public :: lbfgsb_evaluation_error = 53
  integer, parameter, public :: lbfgsb_user_stop = 54

  integer, parameter :: task_new_x = 1
  integer, parameter :: task_start = 2
  integer, parameter :: task_stop = 3
  integer, parameter :: task_fg = 4
  integer, parameter :: task_fg_line_search = 20
  integer, parameter :: task_fg_start = 21
  integer, parameter :: task_x_convergence = 27
  integer, parameter :: task_max_evaluations = 28
  integer, parameter :: task_max_iterations = 29
  integer, parameter :: task_nonfinite = 30
  integer, parameter :: task_invalid_input = 31
  integer, parameter :: task_callback_stop = 32

  type, public :: lbfgsb_control_t
    integer :: memory = 5
    real(dp) :: factr = 1.0e7_dp
    real(dp) :: pgtol = 0.0_dp
    real(dp) :: abstol = 0.0_dp
    real(dp) :: reltol = 1.0e-6_dp
    integer :: max_evaluations = 1000
    integer :: max_iterations = huge(0)
    integer :: trace = 0
    integer :: iprint = -1
    real(dp) :: finite_difference_step = sqrt(epsilon(1.0_dp))
  end type lbfgsb_control_t

  type, public :: lbfgsb_result_t
    real(dp) :: value = huge(1.0_dp)
    real(dp) :: projected_gradient_norm = huge(1.0_dp)
    real(dp), allocatable :: gradient(:)
    integer :: function_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: iterations = 0
    integer :: task = task_invalid_input
    integer :: convergence = lbfgsb_input_error
    integer :: information = 0
    integer :: offending_index = 0
    integer :: icsave = 0
    integer :: lsave(4) = 0
    integer :: isave(44) = 0
    real(dp) :: dsave(29) = 0.0_dp
    logical :: success = .false.
    logical :: initial_point_projected = .false.
    logical :: constrained = .false.
    logical :: all_variables_boxed = .false.
    character(len=:), allocatable :: message
  end type lbfgsb_result_t

  abstract interface
    subroutine objective_gradient_interface(x, f, g, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      class(*), intent(inout), optional :: user_data
    end subroutine objective_gradient_interface

    function objective_interface(x, user_data) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: f
    end function objective_interface

    subroutine monitor_interface(x, f, g, iteration, evaluations, stop, user_data)
      import dp
      real(dp), intent(in) :: x(:), f, g(:)
      integer, intent(in) :: iteration, evaluations
      logical, intent(out) :: stop
      class(*), intent(inout), optional :: user_data
    end subroutine monitor_interface
  end interface

  public :: lbfgsb_minimize
  public :: lbfgsb_minimize_fd
  public :: lbfgsb_task_message

contains

  subroutine lbfgsb_minimize(x, evaluate, result, lower, upper, control, user_data, monitor, function_count_override)
    real(dp), intent(inout) :: x(:)
    procedure(objective_gradient_interface) :: evaluate
    type(lbfgsb_result_t), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(lbfgsb_control_t), intent(in), optional :: control
    class(*), intent(inout), optional :: user_data
    procedure(monitor_interface), optional :: monitor
    integer, intent(in), optional :: function_count_override

    type(lbfgsb_control_t) :: ctrl
    integer :: n, m, nwa, itask, iprint, icsave
    integer :: i, custom_task
    integer, allocatable :: nbd(:), iwa(:)
    integer :: lsave(4), isave(44)
    real(dp) :: f, factr, pgtol
    real(dp), allocatable :: l(:), u(:), g(:), wa(:), last_x(:)
    real(dp) :: dsave(29)
    logical :: stop_requested

    ctrl = lbfgsb_control_t()
    if (present(control)) ctrl = control

    n = size(x)
    call initialize_result(result, n)

    if (.not. valid_control(ctrl, n, result)) return
    if (.not. build_bounds(n, lower, upper, l, u, nbd, result)) return
    if (.not. all(ieee_is_finite(x))) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'initial parameters must be finite')
      return
    end if

    m = ctrl%memory
    nwa = 2*m*n + 11*m*m + 5*n + 8*m
    allocate(wa(nwa), iwa(3*n), g(n), last_x(n))
    wa = 0.0_dp
    iwa = 0
    g = 0.0_dp
    last_x = x
    lsave = 0
    isave = 0
    dsave = 0.0_dp
    icsave = 0
    itask = task_start
    iprint = ctrl%iprint
    call set_core_output(iprint >= 0)
    f = huge(1.0_dp)
    factr = ctrl%factr
    pgtol = ctrl%pgtol
    custom_task = 0

    do
      call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                  itask, iprint, icsave, lsave, isave, dsave)

      select case (itask)
      case (task_fg, task_fg_line_search, task_fg_start)
        if (present(user_data)) then
          call evaluate(x, f, g, user_data)
        else
          call evaluate(x, f, g)
        end if
        if (present(function_count_override)) then
          result%function_evaluations = function_count_override
        else
          result%function_evaluations = result%function_evaluations + 1
        end if
        result%gradient_evaluations = result%gradient_evaluations + 1

        if (.not. ieee_is_finite(f) .or. .not. all(ieee_is_finite(g))) then
          custom_task = task_nonfinite
          exit
        end if

        if (ctrl%trace >= 2) then
          write(output_unit,'(a,i0,a,es16.8,a,es16.8)') &
            'evaluation ', result%function_evaluations, ': f=', f, &
            ', max|g|=', maxval(abs(g))
        end if

      case (task_new_x)
        result%iterations = isave(30)

        if (ctrl%trace >= 1) then
          write(output_unit,'(a,i0,a,es16.8,a,es16.8)') &
            'iteration ', result%iterations, ': f=', f, &
            ', |projected g|=', dsave(13)
        end if

        if (present(monitor)) then
          stop_requested = .false.
          if (present(user_data)) then
            call monitor(x, f, g, result%iterations, &
                         result%function_evaluations, stop_requested, user_data)
          else
            call monitor(x, f, g, result%iterations, &
                         result%function_evaluations, stop_requested)
          end if
          if (stop_requested) then
            custom_task = task_callback_stop
            exit
          end if
        end if

        if (result%function_evaluations > ctrl%max_evaluations) then
          custom_task = task_max_evaluations
          exit
        end if
        if (result%iterations >= ctrl%max_iterations) then
          custom_task = task_max_iterations
          exit
        end if
        if (parameters_converged(last_x, x, ctrl%abstol, ctrl%reltol)) then
          custom_task = task_x_convergence
          exit
        end if
        last_x = x

      case default
        exit
      end select
    end do

    if (custom_task /= 0) itask = custom_task

    result%value = f
    result%gradient = g
    result%task = itask
    result%message = lbfgsb_task_message(itask)
    result%iterations = max(result%iterations, isave(30))
    result%information = isave(35)
    result%offending_index = isave(2)
    result%icsave = icsave
    result%lsave = lsave
    result%isave = isave
    result%dsave = dsave
    result%projected_gradient_norm = dsave(13)
    result%initial_point_projected = lsave(1) /= 0
    result%constrained = lsave(2) /= 0
    result%all_variables_boxed = lsave(3) /= 0
    result%convergence = convergence_code(itask)
    result%success = result%convergence == lbfgsb_success

  end subroutine lbfgsb_minimize

  subroutine lbfgsb_minimize_fd(x, objective, result, lower, upper, control, user_data, monitor)
    real(dp), intent(inout) :: x(:)
    procedure(objective_interface) :: objective
    type(lbfgsb_result_t), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(lbfgsb_control_t), intent(in), optional :: control
    class(*), intent(inout), optional :: user_data
    procedure(monitor_interface), optional :: monitor

    type(lbfgsb_control_t) :: ctrl
    real(dp), allocatable :: lower_work(:), upper_work(:)
    integer :: fd_function_calls

    ctrl = lbfgsb_control_t()
    if (present(control)) ctrl = control
    fd_function_calls = 0

    if (present(lower)) then
      if (size(lower) /= 1 .and. size(lower) /= size(x)) then
        call initialize_result(result, size(x))
        call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                         'lower must have length one or size(x)')
        return
      end if
    end if
    if (present(upper)) then
      if (size(upper) /= 1 .and. size(upper) /= size(x)) then
        call initialize_result(result, size(x))
        call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                         'upper must have length one or size(x)')
        return
      end if
    end if

    allocate(lower_work(size(x)), upper_work(size(x)))
    lower_work = ieee_value(0.0_dp, ieee_negative_inf)
    upper_work = ieee_value(0.0_dp, ieee_positive_inf)
    if (present(lower)) then
      if (size(lower) == 1) then
        lower_work = lower(1)
      else if (size(lower) == size(x)) then
        lower_work = lower
      end if
    end if
    if (present(upper)) then
      if (size(upper) == 1) then
        upper_work = upper(1)
      else if (size(upper) == size(x)) then
        upper_work = upper
      end if
    end if

    if (present(user_data)) then
      if (present(monitor)) then
        call lbfgsb_minimize(x, finite_difference_evaluate, result, &
             lower_work, upper_work, ctrl, user_data, monitor, fd_function_calls)
      else
        call lbfgsb_minimize(x, finite_difference_evaluate, result, &
             lower_work, upper_work, ctrl, user_data, function_count_override=fd_function_calls)
      end if
    else
      if (present(monitor)) then
        call lbfgsb_minimize(x, finite_difference_evaluate, result, &
             lower_work, upper_work, ctrl, monitor=monitor, function_count_override=fd_function_calls)
      else
        call lbfgsb_minimize(x, finite_difference_evaluate, result, &
             lower_work, upper_work, ctrl, function_count_override=fd_function_calls)
      end if
    end if

    result%function_evaluations = fd_function_calls

  contains

    subroutine finite_difference_evaluate(z, f, g, data)
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      class(*), intent(inout), optional :: data
      real(dp), allocatable :: zp(:), zm(:)
      real(dp) :: fp, fm, h, hp, hm
      integer :: j

      allocate(zp(size(z)), zm(size(z)))
      if (present(data)) then
        f = objective(z, data)
      else
        f = objective(z)
      end if
      fd_function_calls = fd_function_calls + 1

      do j = 1, size(z)
        h = ctrl%finite_difference_step * max(1.0_dp, abs(z(j)))
        hp = min(h, upper_work(j) - z(j))
        hm = min(h, z(j) - lower_work(j))

        if (hp > 0.0_dp .and. hm > 0.0_dp) then
          zp = z
          zm = z
          zp(j) = z(j) + hp
          zm(j) = z(j) - hm
          if (present(data)) then
            fp = objective(zp, data)
            fm = objective(zm, data)
          else
            fp = objective(zp)
            fm = objective(zm)
          end if
          fd_function_calls = fd_function_calls + 2
          g(j) = (fp - fm) / (hp + hm)
        else if (hp > 0.0_dp) then
          zp = z
          zp(j) = z(j) + hp
          if (present(data)) then
            fp = objective(zp, data)
          else
            fp = objective(zp)
          end if
          fd_function_calls = fd_function_calls + 1
          g(j) = (fp - f) / hp
        else if (hm > 0.0_dp) then
          zm = z
          zm(j) = z(j) - hm
          if (present(data)) then
            fm = objective(zm, data)
          else
            fm = objective(zm)
          end if
          fd_function_calls = fd_function_calls + 1
          g(j) = (f - fm) / hm
        else
          g(j) = 0.0_dp
        end if
      end do
    end subroutine finite_difference_evaluate

  end subroutine lbfgsb_minimize_fd

  logical function build_bounds(n, lower, upper, l, u, nbd, result) result(ok)
    integer, intent(in) :: n
    real(dp), intent(in), optional :: lower(:), upper(:)
    real(dp), allocatable, intent(out) :: l(:), u(:)
    integer, allocatable, intent(out) :: nbd(:)
    type(lbfgsb_result_t), intent(inout) :: result
    integer :: i

    ok = .false.
    allocate(l(n), u(n), nbd(n))
    l = ieee_value(0.0_dp, ieee_negative_inf)
    u = ieee_value(0.0_dp, ieee_positive_inf)

    if (present(lower)) then
      if (size(lower) == 1) then
        l = lower(1)
      else if (size(lower) == n) then
        l = lower
      else
        call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                         'lower must have length one or size(x)')
        return
      end if
    end if

    if (present(upper)) then
      if (size(upper) == 1) then
        u = upper(1)
      else if (size(upper) == n) then
        u = upper
      else
        call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                         'upper must have length one or size(x)')
        return
      end if
    end if

    do i = 1, n
      if (ieee_is_finite(l(i)) .and. ieee_is_finite(u(i))) then
        nbd(i) = 2
      else if (ieee_is_finite(l(i))) then
        nbd(i) = 1
      else if (ieee_is_finite(u(i))) then
        nbd(i) = 3
      else
        nbd(i) = 0
      end if
      if (nbd(i) == 2 .and. l(i) > u(i)) then
        result%offending_index = i
        call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                         'a lower bound exceeds its upper bound')
        return
      end if
    end do
    ok = .true.
  end function build_bounds

  logical function valid_control(ctrl, n, result) result(ok)
    type(lbfgsb_control_t), intent(in) :: ctrl
    integer, intent(in) :: n
    type(lbfgsb_result_t), intent(inout) :: result

    ok = .false.
    if (n <= 0) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'x must contain at least one parameter')
    else if (ctrl%memory <= 0) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'memory must be positive')
    else if (ctrl%factr < 0.0_dp .or. ctrl%pgtol < 0.0_dp) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'factr and pgtol must be nonnegative')
    else if (ctrl%abstol < 0.0_dp .or. ctrl%reltol < 0.0_dp) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'abstol and reltol must be nonnegative')
    else if (ctrl%max_evaluations <= 0 .or. ctrl%max_iterations <= 0) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'evaluation and iteration limits must be positive')
    else if (ctrl%finite_difference_step <= 0.0_dp) then
      call set_failure(result, task_invalid_input, lbfgsb_input_error, &
                       'finite_difference_step must be positive')
    else
      ok = .true.
    end if
  end function valid_control

  pure logical function parameters_converged(old_x, new_x, abstol, reltol) result(converged)
    real(dp), intent(in) :: old_x(:), new_x(:), abstol, reltol
    if (abstol == 0.0_dp .and. reltol == 0.0_dp) then
      converged = .false.
    else
      converged = all(abs(old_x - new_x) < abstol + reltol*abs(new_x))
    end if
  end function parameters_converged

  subroutine initialize_result(result, n)
    type(lbfgsb_result_t), intent(out) :: result
    integer, intent(in) :: n
    allocate(result%gradient(n))
    result%gradient = 0.0_dp
    result%message = 'not started'
  end subroutine initialize_result

  subroutine set_failure(result, task, convergence, message)
    type(lbfgsb_result_t), intent(inout) :: result
    integer, intent(in) :: task, convergence
    character(len=*), intent(in) :: message
    result%task = task
    result%convergence = convergence
    result%success = .false.
    result%message = message
  end subroutine set_failure

  pure integer function convergence_code(task) result(code)
    integer, intent(in) :: task
    select case (task)
    case (6, 7, 8, task_x_convergence)
      code = lbfgsb_success
    case (task_max_evaluations, task_max_iterations)
      code = lbfgsb_iteration_limit
    case (23:26)
      code = lbfgsb_warning
    case (9:19, task_invalid_input)
      code = lbfgsb_input_error
    case (task_nonfinite)
      code = lbfgsb_evaluation_error
    case (task_callback_stop)
      code = lbfgsb_user_stop
    case default
      code = lbfgsb_warning
    end select
  end function convergence_code

  pure function lbfgsb_task_message(task) result(message)
    integer, intent(in) :: task
    character(len=:), allocatable :: message

    select case (task)
    case (1);  message = 'NEW_X'
    case (2);  message = 'START'
    case (3);  message = 'STOP'
    case (4);  message = 'FG'
    case (5);  message = 'ABNORMAL_TERMINATION_IN_LNSRCH'
    case (6);  message = 'CONVERGENCE'
    case (7);  message = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
    case (8);  message = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
    case (9);  message = 'ERROR: FTOL .LT. ZERO'
    case (10); message = 'ERROR: GTOL .LT. ZERO'
    case (11); message = 'ERROR: INITIAL G .GE. ZERO'
    case (12); message = 'ERROR: INVALID NBD OR MEMORY'
    case (13); message = 'ERROR: N .LE. 0'
    case (14); message = 'ERROR: NO FEASIBLE SOLUTION'
    case (15); message = 'ERROR: STP .GT. STPMAX'
    case (16); message = 'ERROR: STP .LT. STPMIN'
    case (17); message = 'ERROR: STPMAX .LT. STPMIN'
    case (18); message = 'ERROR: STPMIN .LT. ZERO'
    case (19); message = 'ERROR: XTOL .LT. ZERO'
    case (20); message = 'FG_LNSRCH'
    case (21); message = 'FG_START'
    case (22); message = 'RESTART_FROM_LNSRCH'
    case (23); message = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
    case (24); message = 'WARNING: STP = STPMAX'
    case (25); message = 'WARNING: STP = STPMIN'
    case (26); message = 'WARNING: XTOL TEST SATISFIED'
    case (task_x_convergence)
      message = 'CONVERGENCE: PARAMETER DIFFERENCES BELOW XTOL'
    case (task_max_evaluations)
      message = 'MAXIMUM NUMBER OF FUNCTION EVALUATIONS REACHED'
    case (task_max_iterations)
      message = 'MAXIMUM NUMBER OF ITERATIONS REACHED'
    case (task_nonfinite)
      message = 'OBJECTIVE OR GRADIENT RETURNED A NON-FINITE VALUE'
    case (task_invalid_input)
      message = 'INVALID INPUT'
    case (task_callback_stop)
      message = 'STOPPED BY MONITOR CALLBACK'
    case default
      message = 'UNKNOWN L-BFGS-B TASK'
    end select
  end function lbfgsb_task_message

end module lbfgsb3_mod
