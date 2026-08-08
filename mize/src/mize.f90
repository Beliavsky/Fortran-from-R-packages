! Modern Fortran computational port of mize.
! Copyright (c) 2016 James Melville
! Copyright (c) 2026 Fortran port contributors
! SPDX-License-Identifier: BSD-2-Clause

module mize_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)

  integer, parameter, public :: mize_success = 0
  integer, parameter, public :: mize_max_iterations = 1
  integer, parameter, public :: mize_max_evaluations = 2
  integer, parameter, public :: mize_line_search_failure = 3
  integer, parameter, public :: mize_invalid_input = 10
  integer, parameter, public :: mize_nonfinite = 11
  integer, parameter, public :: mize_user_stop = 12

  type, public :: mize_control_t
    character(len=24) :: method = 'l-bfgs'
    character(len=24) :: line_search = 'more-thuente'
    character(len=16) :: cg_update = 'pr+'
    character(len=16) :: preconditioner = ''
    character(len=16) :: tn_init = 'zero'
    character(len=16) :: tn_exit = 'curvature'
    character(len=16) :: mom_type = 'classical'
    character(len=16) :: mom_schedule = 'constant'
    character(len=16) :: restart = ''
    integer :: memory = 5
    integer :: max_iterations = 100
    integer :: max_evaluations = 100000
    integer :: max_line_search = 40
    integer :: tn_max_iterations = 40
    integer :: hessian_every = 0
    integer :: restart_wait = 10
    real(dp) :: step0 = 1.0_dp
    real(dp) :: c1 = 1.0e-4_dp
    real(dp) :: c2 = 0.9_dp
    real(dp) :: min_step = 1.0e-20_dp
    real(dp) :: max_step = 1.0e20_dp
    real(dp) :: step_up = 1.1_dp
    real(dp) :: step_down = 0.5_dp
    real(dp) :: dbd_weight = 0.1_dp
    real(dp) :: mom_init = 0.9_dp
    real(dp) :: mom_final = 0.9_dp
    integer :: mom_switch_iter = 50
    real(dp) :: nest_q = 0.0_dp
    integer :: nest_burn_in = 0
    real(dp) :: abs_tol = sqrt(epsilon(1.0_dp))
    real(dp) :: rel_tol = sqrt(epsilon(1.0_dp))
    real(dp) :: grad_tol = sqrt(epsilon(1.0_dp))
    real(dp) :: ginf_tol = sqrt(epsilon(1.0_dp))
    real(dp) :: step_tol = sqrt(epsilon(1.0_dp))
    logical :: scale_hess = .true.
    logical :: norm_direction = .false.
    logical :: strong_curvature = .true.
    logical :: approx_armijo = .false.
    logical :: step_up_add = .false.
    logical :: mom_linear_weight = .false.
    logical :: use_init_momentum = .false.
    logical :: store_progress = .false.
  end type mize_control_t

  type, public :: mize_result_t
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: best_par(:)
    real(dp) :: value = huge(1.0_dp)
    real(dp) :: best_value = huge(1.0_dp)
    real(dp) :: gradient_norm = huge(1.0_dp)
    real(dp) :: gradient_inf_norm = huge(1.0_dp)
    integer :: iterations = 0
    integer :: function_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: status = mize_invalid_input
    logical :: converged = .false.
    character(len=:), allocatable :: message
    real(dp), allocatable :: progress_value(:)
    real(dp), allocatable :: progress_gradient_norm(:)
    real(dp), allocatable :: progress_step_norm(:)
  end type mize_result_t

  type, public :: mize_state_t
    type(mize_control_t) :: control
    logical :: initialized = .false.
    logical :: terminated = .false.
    logical :: have_previous = .false.
    logical :: have_fixed_hessian = .false.
    integer :: n = 0
    integer :: iteration = 0
    integer :: function_evaluations = 0
    integer :: gradient_evaluations = 0
    integer :: status = mize_invalid_input
    integer :: last_restart = -1000000
    character(len=:), allocatable :: message
    real(dp) :: value = huge(1.0_dp)
    real(dp) :: previous_value = huge(1.0_dp)
    real(dp) :: best_value = huge(1.0_dp)
    real(dp) :: last_step = 0.0_dp
    real(dp) :: last_step_norm = huge(1.0_dp)
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: gradient(:)
    real(dp), allocatable :: previous_gradient(:)
    real(dp), allocatable :: previous_direction(:)
    real(dp), allocatable :: previous_update(:)
    real(dp), allocatable :: best_par(:)
    real(dp), allocatable :: inverse_hessian(:, :)
    real(dp), allocatable :: fixed_hessian(:, :)
    real(dp), allocatable :: s_history(:, :)
    real(dp), allocatable :: y_history(:, :)
    real(dp), allocatable :: rho_history(:)
    integer :: history_count = 0
    real(dp), allocatable :: velocity(:)
    real(dp), allocatable :: dbd_average(:)
    real(dp), allocatable :: dbd_scale(:)
    real(dp), allocatable :: progress_value(:)
    real(dp), allocatable :: progress_gradient_norm(:)
    real(dp), allocatable :: progress_step_norm(:)
  end type mize_state_t

  type, public :: gradient_check_t
    real(dp), allocatable :: analytic(:)
    real(dp), allocatable :: finite_difference(:)
    real(dp), allocatable :: absolute_error(:)
    real(dp), allocatable :: relative_error(:)
    real(dp) :: maximum_absolute_error = 0.0_dp
    real(dp) :: maximum_relative_error = 0.0_dp
    integer :: worst_coordinate = 0
    logical :: passed = .false.
  end type gradient_check_t

  abstract interface
    subroutine mize_fg_interface(x, f, g, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      class(*), intent(inout), optional :: user_data
    end subroutine mize_fg_interface

    subroutine mize_hessian_interface(x, h, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: h(:, :)
      class(*), intent(inout), optional :: user_data
    end subroutine mize_hessian_interface

    subroutine mize_hvp_interface(x, v, hv, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: hv(:)
      class(*), intent(inout), optional :: user_data
    end subroutine mize_hvp_interface

    subroutine mize_monitor_interface(x, f, g, iteration, function_evaluations, &
                                      gradient_evaluations, stop, user_data)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f
      real(dp), intent(in) :: g(:)
      integer, intent(in) :: iteration
      integer, intent(in) :: function_evaluations
      integer, intent(in) :: gradient_evaluations
      logical, intent(out) :: stop
      class(*), intent(inout), optional :: user_data
    end subroutine mize_monitor_interface

    subroutine mize_momentum_interface(iteration, max_iterations, momentum, user_data)
      import dp
      integer, intent(in) :: iteration
      integer, intent(in) :: max_iterations
      real(dp), intent(out) :: momentum
      class(*), intent(inout), optional :: user_data
    end subroutine mize_momentum_interface
  end interface

  public :: mize_minimize
  public :: mize_init
  public :: mize_step
  public :: mize_state_result
  public :: check_mize_convergence
  public :: check_mize_gradient
  public :: mize_status_message

contains

  subroutine mize_minimize(par, fg, result, control, hessian, hvp, monitor, &
                           momentum_schedule, user_data)
    real(dp), intent(inout) :: par(:)
    procedure(mize_fg_interface) :: fg
    type(mize_result_t), intent(out) :: result
    type(mize_control_t), intent(in), optional :: control
    procedure(mize_hessian_interface), optional :: hessian
    procedure(mize_hvp_interface), optional :: hvp
    procedure(mize_monitor_interface), optional :: monitor
    procedure(mize_momentum_interface), optional :: momentum_schedule
    class(*), intent(inout), optional :: user_data

    type(mize_state_t) :: state

    call mize_init(state, par, fg, control, user_data)
    if (.not. state%initialized) then
      call mize_state_result(state, result)
      return
    end if

    do while (.not. state%terminated)
      call mize_step(state, fg, hessian, hvp, monitor, momentum_schedule, user_data)
    end do

    par = state%par
    call mize_state_result(state, result)
  end subroutine mize_minimize

  subroutine mize_init(state, par, fg, control, user_data)
    type(mize_state_t), intent(out) :: state
    real(dp), intent(in) :: par(:)
    procedure(mize_fg_interface) :: fg
    type(mize_control_t), intent(in), optional :: control
    class(*), intent(inout), optional :: user_data

    integer :: n, i
    logical :: ok

    state%control = mize_control_t()
    if (present(control)) state%control = control

    n = size(par)
    if (n < 1) then
      call fail_state(state, mize_invalid_input, 'parameter vector must not be empty')
      return
    end if
    if (.not. all(ieee_is_finite(par))) then
      call fail_state(state, mize_invalid_input, 'initial parameters must be finite')
      return
    end if
    if (.not. valid_control(state%control)) then
      call fail_state(state, mize_invalid_input, 'invalid optimizer control values')
      return
    end if

    state%n = n
    allocate(state%par(n), state%gradient(n), state%previous_gradient(n))
    allocate(state%previous_direction(n), state%previous_update(n), state%best_par(n))
    allocate(state%inverse_hessian(n, n), state%fixed_hessian(n, n))
    allocate(state%s_history(n, state%control%memory))
    allocate(state%y_history(n, state%control%memory))
    allocate(state%rho_history(state%control%memory))
    allocate(state%velocity(n), state%dbd_average(n), state%dbd_scale(n))

    state%par = par
    state%previous_gradient = 0.0_dp
    state%previous_direction = 0.0_dp
    state%previous_update = 0.0_dp
    state%inverse_hessian = 0.0_dp
    state%fixed_hessian = 0.0_dp
    do i = 1, n
      state%inverse_hessian(i, i) = 1.0_dp
    end do
    state%s_history = 0.0_dp
    state%y_history = 0.0_dp
    state%rho_history = 0.0_dp
    state%velocity = 0.0_dp
    state%dbd_average = 0.0_dp
    state%dbd_scale = state%control%step0

    if (state%control%store_progress) then
      allocate(state%progress_value(0:state%control%max_iterations))
      allocate(state%progress_gradient_norm(0:state%control%max_iterations))
      allocate(state%progress_step_norm(0:state%control%max_iterations))
      state%progress_value = huge(1.0_dp)
      state%progress_gradient_norm = huge(1.0_dp)
      state%progress_step_norm = huge(1.0_dp)
    end if

    call evaluate_fg(fg, state%par, state%value, state%gradient, user_data, ok)
    state%function_evaluations = 1
    state%gradient_evaluations = 1
    if (.not. ok) then
      call fail_state(state, mize_nonfinite, 'initial objective or gradient is non-finite')
      return
    end if

    state%best_value = state%value
    state%best_par = state%par
    state%previous_value = state%value
    state%initialized = .true.
    state%status = mize_success
    state%message = 'initialized'

    if (state%control%store_progress) then
      state%progress_value(0) = state%value
      state%progress_gradient_norm(0) = vector_norm(state%gradient)
      state%progress_step_norm(0) = 0.0_dp
    end if

    if (gradient_converged(state)) then
      state%terminated = .true.
      state%status = mize_success
      state%message = 'initial gradient satisfies convergence tolerance'
    end if
  end subroutine mize_init

  subroutine mize_step(state, fg, hessian, hvp, monitor, momentum_schedule, user_data)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    procedure(mize_hessian_interface), optional :: hessian
    procedure(mize_hvp_interface), optional :: hvp
    procedure(mize_monitor_interface), optional :: monitor
    procedure(mize_momentum_interface), optional :: momentum_schedule
    class(*), intent(inout), optional :: user_data

    real(dp), allocatable :: x_old(:), g_old(:), direction(:), update(:)
    real(dp), allocatable :: x_new(:), g_new(:), s(:), y(:)
    real(dp) :: f_old, f_new, alpha
    integer :: ls_status
    logical :: ok, stop
    character(len=:), allocatable :: method

    if (.not. state%initialized .or. state%terminated) return
    if (state%iteration >= state%control%max_iterations) then
      call terminate_state(state, mize_max_iterations, 'maximum iterations reached')
      return
    end if
    if (state%function_evaluations >= state%control%max_evaluations) then
      call terminate_state(state, mize_max_evaluations, 'maximum evaluations reached')
      return
    end if

    allocate(x_old(state%n), g_old(state%n), direction(state%n), update(state%n))
    allocate(x_new(state%n), g_new(state%n), s(state%n), y(state%n))
    x_old = state%par
    g_old = state%gradient
    f_old = state%value
    method = canonical_name(state%control%method)

    if (method == 'dbd') then
      call dbd_step(state, fg, x_new, f_new, g_new, user_data, ok)
      if (.not. ok) then
        call terminate_state(state, mize_nonfinite, 'delta-bar-delta produced non-finite values')
        return
      end if
      update = x_new - x_old
      alpha = 1.0_dp
    else
      call calculate_direction(state, fg, direction, hessian, hvp, momentum_schedule, user_data, ok)
      if (.not. ok) then
        call terminate_state(state, mize_nonfinite, 'failed to calculate a finite search direction')
        return
      end if
      if (state%control%norm_direction) call normalize_vector(direction)
      if (dot_product(g_old, direction) >= 0.0_dp) direction = -g_old

      call perform_line_search(state, fg, direction, x_new, f_new, g_new, alpha, &
                               ls_status, user_data)
      if (ls_status /= mize_success) then
        call terminate_state(state, ls_status, mize_status_message(ls_status))
        return
      end if
      update = x_new - x_old
    end if

    s = update
    y = g_new - g_old
    call update_method_state(state, method, s, y)
    call apply_adaptive_restart(state, method, f_old, f_new, g_new, update)

    state%previous_value = f_old
    state%previous_gradient = g_old
    state%previous_direction = direction
    state%previous_update = update
    state%have_previous = .true.
    state%par = x_new
    state%value = f_new
    state%gradient = g_new
    state%last_step = alpha
    state%last_step_norm = vector_norm(update)
    state%iteration = state%iteration + 1

    if (f_new < state%best_value) then
      state%best_value = f_new
      state%best_par = x_new
    end if

    if (state%control%store_progress) then
      state%progress_value(state%iteration) = state%value
      state%progress_gradient_norm(state%iteration) = vector_norm(state%gradient)
      state%progress_step_norm(state%iteration) = state%last_step_norm
    end if

    call check_mize_convergence(state)

    if (present(monitor) .and. .not. state%terminated) then
      stop = .false.
      call invoke_monitor(monitor, state%par, state%value, state%gradient, &
                          state%iteration, state%function_evaluations, &
                          state%gradient_evaluations, stop, user_data)
      if (stop) call terminate_state(state, mize_user_stop, 'stopped by monitor callback')
    end if
  end subroutine mize_step

  subroutine calculate_direction(state, fg, direction, hessian, hvp, momentum_schedule, user_data, ok)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(out) :: direction(:)
    procedure(mize_hessian_interface), optional :: hessian
    procedure(mize_hvp_interface), optional :: hvp
    procedure(mize_momentum_interface), optional :: momentum_schedule
    class(*), intent(inout), optional :: user_data
    logical, intent(out) :: ok

    character(len=:), allocatable :: method
    real(dp), allocatable :: h(:, :), work(:), lookahead(:), look_g(:)
    real(dp) :: look_f, mu, denom
    logical :: eval_ok, solve_ok

    method = canonical_name(state%control%method)
    ok = .true.

    select case (method)
    case ('sd')
      direction = -state%gradient

    case ('bfgs', 'sr1')
      direction = -matmul(state%inverse_hessian, state%gradient)

    case ('l-bfgs')
      call lbfgs_apply(state, state%gradient, direction)
      direction = -direction

    case ('cg')
      call cg_direction_value(state, direction)

    case ('newton')
      allocate(h(state%n, state%n))
      call obtain_hessian(state, fg, h, hessian, user_data, eval_ok)
      if (.not. eval_ok) then
        ok = .false.
        return
      end if
      allocate(work(state%n))
      work = -state%gradient
      call solve_linear_system(h, work, direction, solve_ok)
      if (.not. solve_ok) direction = -state%gradient

    case ('phess')
      if (.not. state%have_fixed_hessian) then
        call obtain_hessian(state, fg, state%fixed_hessian, hessian, user_data, eval_ok)
        if (.not. eval_ok) then
          ok = .false.
          return
        end if
        state%have_fixed_hessian = .true.
      else if (state%control%hessian_every > 0) then
        if (mod(state%iteration, state%control%hessian_every) == 0) then
          call obtain_hessian(state, fg, state%fixed_hessian, hessian, user_data, eval_ok)
          if (.not. eval_ok) then
            ok = .false.
            return
          end if
        end if
      end if
      allocate(work(state%n))
      work = -state%gradient
      call solve_linear_system(state%fixed_hessian, work, direction, solve_ok)
      if (.not. solve_ok) direction = -state%gradient

    case ('tn')
      call truncated_newton_direction(state, fg, direction, hessian, hvp, user_data, eval_ok)
      if (.not. eval_ok) then
        ok = .false.
        return
      end if

    case ('momentum')
      mu = momentum_value(state, momentum_schedule, user_data)
      if (state%iteration == 0 .and. .not. state%control%use_init_momentum) mu = 0.0_dp
      if (canonical_name(state%control%mom_type) == 'nesterov') then
        allocate(lookahead(state%n), look_g(state%n))
        lookahead = state%par + mu * state%velocity
        call evaluate_fg(fg, lookahead, look_f, look_g, user_data, eval_ok)
        state%function_evaluations = state%function_evaluations + 1
        state%gradient_evaluations = state%gradient_evaluations + 1
        if (.not. eval_ok) then
          ok = .false.
          return
        end if
        if (state%control%mom_linear_weight) then
          direction = -(1.0_dp - mu) * look_g + mu * state%velocity / &
                      max(state%control%step0, tiny(1.0_dp))
        else
          direction = -look_g + mu * state%velocity / &
                      max(state%control%step0, tiny(1.0_dp))
        end if
      else
        if (state%control%mom_linear_weight) then
          direction = -(1.0_dp - mu) * state%gradient + mu * state%velocity / &
                      max(state%control%step0, tiny(1.0_dp))
        else
          direction = -state%gradient + mu * state%velocity / &
                      max(state%control%step0, tiny(1.0_dp))
        end if
      end if

    case ('nag')
      mu = nesterov_momentum(state%iteration + 1, state%control%nest_q, &
                             state%control%nest_burn_in)
      if (state%iteration == 0 .and. .not. state%control%use_init_momentum) mu = 0.0_dp
      allocate(lookahead(state%n), look_g(state%n))
      lookahead = state%par + mu * state%velocity
      call evaluate_fg(fg, lookahead, look_f, look_g, user_data, eval_ok)
      state%function_evaluations = state%function_evaluations + 1
      state%gradient_evaluations = state%gradient_evaluations + 1
      if (.not. eval_ok) then
        ok = .false.
        return
      end if
      denom = max(state%control%step0, tiny(1.0_dp))
      direction = -look_g + mu * state%velocity / denom

    case default
      ok = .false.
      direction = 0.0_dp
    end select

    if (.not. all(ieee_is_finite(direction))) ok = .false.
  end subroutine calculate_direction

  subroutine perform_line_search(state, fg, direction, x_new, f_new, g_new, alpha, &
                                 status, user_data)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: direction(:)
    real(dp), intent(out) :: x_new(:)
    real(dp), intent(out) :: f_new
    real(dp), intent(out) :: g_new(:)
    real(dp), intent(out) :: alpha
    integer, intent(out) :: status
    class(*), intent(inout), optional :: user_data

    character(len=:), allocatable :: ls
    real(dp) :: d0

    ls = canonical_name(state%control%line_search)
    d0 = dot_product(state%gradient, direction)
    if (.not. ieee_is_finite(d0) .or. d0 >= 0.0_dp) then
      status = mize_line_search_failure
      return
    end if

    select case (ls)
    case ('constant')
      alpha = max(state%control%min_step, min(state%control%max_step, &
                  state%control%step0))
      call evaluate_trial(state, fg, direction, alpha, x_new, f_new, g_new, status, user_data)

    case ('backtracking', 'armijo')
      call armijo_backtracking(state, fg, direction, .true., x_new, f_new, g_new, &
                               alpha, status, user_data)

    case ('bold', 'bold-driver', 'bold driver')
      call armijo_backtracking(state, fg, direction, .false., x_new, f_new, g_new, &
                               alpha, status, user_data)

    case default
      call wolfe_line_search(state, fg, direction, x_new, f_new, g_new, alpha, &
                             status, user_data)
    end select
  end subroutine perform_line_search

  subroutine armijo_backtracking(state, fg, direction, use_armijo, x_new, f_new, &
                                 g_new, alpha, status, user_data)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: direction(:)
    logical, intent(in) :: use_armijo
    real(dp), intent(out) :: x_new(:), f_new, g_new(:), alpha
    integer, intent(out) :: status
    class(*), intent(inout), optional :: user_data

    integer :: i
    real(dp) :: d0, rhs

    d0 = dot_product(state%gradient, direction)
    if (canonical_name(state%control%line_search) == 'bold' .or. &
        canonical_name(state%control%line_search) == 'bold-driver' .or. &
        canonical_name(state%control%line_search) == 'bold driver') then
      if (state%last_step > 0.0_dp) then
        alpha = min(state%control%max_step, state%last_step * state%control%step_up)
      else
        alpha = state%control%step0
      end if
    else
      alpha = initial_step(state)
    end if

    do i = 1, state%control%max_line_search
      call evaluate_trial(state, fg, direction, alpha, x_new, f_new, g_new, status, user_data)
      if (status == mize_success) then
        if (use_armijo) then
          rhs = state%value + state%control%c1 * alpha * d0
          if (f_new <= rhs) return
        else
          if (f_new < state%value) return
        end if
      end if
      alpha = alpha * state%control%step_down
      if (alpha < state%control%min_step) exit
      if (state%function_evaluations >= state%control%max_evaluations) exit
    end do
    status = mize_line_search_failure
  end subroutine armijo_backtracking

  subroutine wolfe_line_search(state, fg, direction, x_new, f_new, g_new, alpha, &
                               status, user_data)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: direction(:)
    real(dp), intent(out) :: x_new(:), f_new, g_new(:), alpha
    integer, intent(out) :: status
    class(*), intent(inout), optional :: user_data

    real(dp), allocatable :: x_prev(:), g_prev(:)
    real(dp) :: alpha_prev, f_prev, d_prev, d0, dnew
    integer :: i
    logical :: armijo_ok, curvature_ok

    allocate(x_prev(state%n), g_prev(state%n))
    alpha_prev = 0.0_dp
    f_prev = state%value
    d0 = dot_product(state%gradient, direction)
    d_prev = d0
    x_prev = state%par
    g_prev = state%gradient
    alpha = initial_step(state)

    do i = 1, state%control%max_line_search
      call evaluate_trial(state, fg, direction, alpha, x_new, f_new, g_new, status, user_data)
      if (status /= mize_success) then
        alpha = 0.5_dp * (alpha_prev + alpha)
        if (alpha <= state%control%min_step) exit
        cycle
      end if
      dnew = dot_product(g_new, direction)
      armijo_ok = line_armijo_ok(state, f_new, dnew, alpha, d0)
      if ((.not. armijo_ok) .or. (i > 1 .and. f_new >= f_prev)) then
        block
          real(dp) :: alo_call, ahi_call, flo_call, dlo_call
          alo_call = alpha_prev
          ahi_call = alpha
          flo_call = f_prev
          dlo_call = d_prev
          call zoom_line_search(state, fg, direction, alo_call, ahi_call, flo_call, &
                                dlo_call, x_new, f_new, g_new, alpha, status, user_data)
        end block
        return
      end if

      if (state%control%strong_curvature) then
        curvature_ok = abs(dnew) <= -state%control%c2 * d0
      else
        curvature_ok = dnew >= state%control%c2 * d0
      end if
      if (curvature_ok) then
        status = mize_success
        return
      end if
      if (dnew >= 0.0_dp) then
        block
          real(dp) :: alo_call, ahi_call, flo_call, dlo_call
          alo_call = alpha
          ahi_call = alpha_prev
          flo_call = f_new
          dlo_call = dnew
          call zoom_line_search(state, fg, direction, alo_call, ahi_call, flo_call, &
                                dlo_call, x_new, f_new, g_new, alpha, status, user_data)
        end block
        return
      end if

      alpha_prev = alpha
      f_prev = f_new
      d_prev = dnew
      x_prev = x_new
      g_prev = g_new
      alpha = min(state%control%max_step, 2.0_dp * alpha)
      if (alpha <= alpha_prev + state%control%min_step) exit
      if (state%function_evaluations >= state%control%max_evaluations) exit
    end do

    call armijo_backtracking(state, fg, direction, .true., x_new, f_new, g_new, &
                             alpha, status, user_data)
  end subroutine wolfe_line_search

  subroutine zoom_line_search(state, fg, direction, alo_in, ahi_in, flo_in, dlo_in, &
                              x_new, f_new, g_new, alpha, status, user_data)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: direction(:)
    real(dp), intent(in) :: alo_in, ahi_in, flo_in, dlo_in
    real(dp), intent(out) :: x_new(:), f_new, g_new(:), alpha
    integer, intent(out) :: status
    class(*), intent(inout), optional :: user_data

    real(dp) :: alo, ahi, flo, fhi, dlo, dhi, d0, da
    real(dp), allocatable :: xt(:), gt(:)
    integer :: i, trial_status
    logical :: curvature_ok

    allocate(xt(state%n), gt(state%n))
    alo = alo_in
    ahi = ahi_in
    flo = flo_in
    dlo = dlo_in
    d0 = dot_product(state%gradient, direction)

    if (abs(ahi) <= tiny(1.0_dp)) then
      fhi = state%value
      dhi = d0
    else
      call evaluate_trial(state, fg, direction, ahi, xt, fhi, gt, trial_status, user_data)
      if (trial_status /= mize_success) then
        fhi = huge(1.0_dp)
        dhi = huge(1.0_dp)
      else
        dhi = dot_product(gt, direction)
      end if
    end if

    do i = 1, state%control%max_line_search
      alpha = safeguarded_interpolation(alo, flo, dlo, ahi, fhi, dhi)
      call evaluate_trial(state, fg, direction, alpha, x_new, f_new, g_new, status, user_data)
      if (status /= mize_success) then
        ahi = alpha
        fhi = huge(1.0_dp)
        dhi = huge(1.0_dp)
        cycle
      end if
      da = dot_product(g_new, direction)
      if ((.not. line_armijo_ok(state, f_new, da, alpha, d0)) .or. f_new >= flo) then
        ahi = alpha
        fhi = f_new
        dhi = da
      else
        if (state%control%strong_curvature) then
          curvature_ok = abs(da) <= -state%control%c2 * d0
        else
          curvature_ok = da >= state%control%c2 * d0
        end if
        if (curvature_ok) then
          status = mize_success
          return
        end if
        if (da * (ahi - alo) >= 0.0_dp) then
          ahi = alo
          fhi = flo
          dhi = dlo
        end if
        alo = alpha
        flo = f_new
        dlo = da
      end if
      if (abs(ahi - alo) <= state%control%min_step * max(1.0_dp, abs(alo))) then
        status = mize_success
        return
      end if
      if (state%function_evaluations >= state%control%max_evaluations) exit
    end do
    status = mize_line_search_failure
  end subroutine zoom_line_search

  logical function line_armijo_ok(state, f, d, alpha, d0) result(ok)
    type(mize_state_t), intent(in) :: state
    real(dp), intent(in) :: f, d, alpha, d0
    real(dp) :: eps_k

    ok = f <= state%value + state%control%c1 * alpha * d0
    if (.not. ok .and. state%control%approx_armijo) then
      eps_k = 1.0e-6_dp * abs(state%value)
      ok = f <= state%value + eps_k .and. &
           d <= (2.0_dp * state%control%c1 - 1.0_dp) * d0
    end if
  end function line_armijo_ok

  real(dp) function safeguarded_interpolation(a, fa, da, b, fb, db) result(x)
    real(dp), intent(in) :: a, fa, da, b, fb, db
    real(dp) :: lo, hi, d1, rad, denom, xc

    lo = min(a, b)
    hi = max(a, b)
    x = 0.5_dp * (a + b)
    if (.not. all(ieee_is_finite([fa, da, fb, db]))) return
    d1 = da + db - 3.0_dp * (fa - fb) / (a - b)
    rad = d1 * d1 - da * db
    if (rad < 0.0_dp .or. .not. ieee_is_finite(rad)) return
    rad = sqrt(rad)
    denom = db - da + 2.0_dp * rad
    if (abs(denom) <= tiny(1.0_dp)) return
    if (b > a) then
      xc = b - (b - a) * (db + rad - d1) / denom
    else
      xc = b - (b - a) * (db - rad - d1) / (db - da - 2.0_dp * rad)
    end if
    if (ieee_is_finite(xc) .and. xc > lo + 0.1_dp * (hi - lo) .and. &
        xc < hi - 0.1_dp * (hi - lo)) x = xc
  end function safeguarded_interpolation

  subroutine evaluate_trial(state, fg, direction, alpha, x, f, g, status, user_data)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: direction(:), alpha
    real(dp), intent(out) :: x(:), f, g(:)
    integer, intent(out) :: status
    class(*), intent(inout), optional :: user_data
    logical :: ok

    x = state%par + alpha * direction
    call evaluate_fg(fg, x, f, g, user_data, ok)
    state%function_evaluations = state%function_evaluations + 1
    state%gradient_evaluations = state%gradient_evaluations + 1
    if (ok) then
      status = mize_success
    else
      status = mize_nonfinite
    end if
  end subroutine evaluate_trial

  subroutine update_method_state(state, method, s, y)
    type(mize_state_t), intent(inout) :: state
    character(len=*), intent(in) :: method
    real(dp), intent(in) :: s(:), y(:)
    real(dp) :: ys, yy, rho, scale
    real(dp), allocatable :: hy(:), v(:), hnew(:, :)
    integer :: i, j

    select case (method)
    case ('bfgs')
      ys = dot_product(y, s)
      if (curvature_ok_update(y, s)) then
        if (state%iteration == 0 .and. state%control%scale_hess) then
          yy = dot_product(y, y)
          if (yy > tiny(1.0_dp)) state%inverse_hessian = (ys / yy) * state%inverse_hessian
        end if
        rho = 1.0_dp / ys
        allocate(hnew(state%n, state%n), hy(state%n))
        hy = matmul(state%inverse_hessian, y)
        scale = (1.0_dp + rho * dot_product(y, hy)) * rho
        do j = 1, state%n
          do i = 1, state%n
            hnew(i, j) = state%inverse_hessian(i, j) + scale * s(i) * s(j) - &
                         rho * (hy(i) * s(j) + s(i) * hy(j))
          end do
        end do
        state%inverse_hessian = hnew
      end if

    case ('sr1')
      allocate(v(state%n))
      v = s - matmul(state%inverse_hessian, y)
      ys = dot_product(v, y)
      if (abs(ys) > 1.0e-8_dp * vector_norm(v) * vector_norm(y)) then
        do j = 1, state%n
          do i = 1, state%n
            state%inverse_hessian(i, j) = state%inverse_hessian(i, j) + &
                                          v(i) * v(j) / ys
          end do
        end do
      end if
      if (dot_product(state%gradient, -matmul(state%inverse_hessian, state%gradient)) >= 0.0_dp) then
        state%inverse_hessian = 0.0_dp
        do i = 1, state%n
          state%inverse_hessian(i, i) = 1.0_dp
        end do
      end if

    case ('l-bfgs', 'cg', 'tn')
      call add_lbfgs_history(state, s, y)

    case ('momentum', 'nag')
      state%velocity = s

    case default
      continue
    end select
  end subroutine update_method_state

  subroutine add_lbfgs_history(state, s, y)
    type(mize_state_t), intent(inout) :: state
    real(dp), intent(in) :: s(:), y(:)
    real(dp) :: ys
    integer :: m

    if (.not. curvature_ok_update(y, s)) return
    ys = dot_product(y, s)
    m = state%control%memory
    if (state%history_count < m) then
      state%history_count = state%history_count + 1
      state%s_history(:, state%history_count) = s
      state%y_history(:, state%history_count) = y
      state%rho_history(state%history_count) = 1.0_dp / ys
    else
      state%s_history(:, 1:m-1) = state%s_history(:, 2:m)
      state%y_history(:, 1:m-1) = state%y_history(:, 2:m)
      state%rho_history(1:m-1) = state%rho_history(2:m)
      state%s_history(:, m) = s
      state%y_history(:, m) = y
      state%rho_history(m) = 1.0_dp / ys
    end if
  end subroutine add_lbfgs_history

  subroutine lbfgs_apply(state, q_in, z)
    type(mize_state_t), intent(in) :: state
    real(dp), intent(in) :: q_in(:)
    real(dp), intent(out) :: z(:)
    real(dp), allocatable :: q(:), alpha(:)
    real(dp) :: beta, gamma, ys, yy
    integer :: i

    allocate(q(state%n), alpha(max(1, state%history_count)))
    q = q_in
    if (state%history_count == 0) then
      z = q
      return
    end if

    do i = state%history_count, 1, -1
      alpha(i) = state%rho_history(i) * dot_product(state%s_history(:, i), q)
      q = q - alpha(i) * state%y_history(:, i)
    end do

    gamma = 1.0_dp
    if (state%control%scale_hess) then
      ys = dot_product(state%s_history(:, state%history_count), &
                       state%y_history(:, state%history_count))
      yy = dot_product(state%y_history(:, state%history_count), &
                       state%y_history(:, state%history_count))
      if (yy > tiny(1.0_dp)) gamma = ys / yy
    end if
    z = gamma * q

    do i = 1, state%history_count
      beta = state%rho_history(i) * dot_product(state%y_history(:, i), z)
      z = z + state%s_history(:, i) * (alpha(i) - beta)
    end do
  end subroutine lbfgs_apply

  subroutine cg_direction_value(state, direction)
    type(mize_state_t), intent(in) :: state
    real(dp), intent(out) :: direction(:)
    real(dp), allocatable :: y(:), w(:), wold(:), v(:)
    real(dp) :: beta, denom, eta_k, eps
    character(len=:), allocatable :: update

    if (.not. state%have_previous) then
      direction = -state%gradient
      return
    end if

    allocate(y(state%n), w(state%n), wold(state%n), v(state%n))
    y = state%gradient - state%previous_gradient
    w = state%gradient
    wold = state%previous_gradient
    if (canonical_name(state%control%preconditioner) == 'l-bfgs' .and. &
        state%history_count > 0) then
      call lbfgs_apply(state, state%gradient, w)
      call lbfgs_apply(state, state%previous_gradient, wold)
    end if

    eps = epsilon(1.0_dp)
    update = canonical_name(state%control%cg_update)
    select case (update)
    case ('fr')
      beta = dot_product(state%gradient, w) / &
             (dot_product(state%previous_gradient, wold) + eps)
    case ('cd')
      beta = -dot_product(state%gradient, w) / &
             (dot_product(state%previous_direction, state%previous_gradient) + eps)
    case ('dy')
      beta = dot_product(state%gradient, w) / &
             (dot_product(state%previous_direction, y) + eps)
    case ('hs', 'hs+')
      beta = dot_product(y, w) / (dot_product(state%previous_direction, y) + eps)
      if (update == 'hs+') beta = max(0.0_dp, beta)
    case ('pr', 'pr+')
      beta = dot_product(w, y) / &
             (dot_product(state%previous_gradient, wold) + eps)
      if (update == 'pr+') beta = max(0.0_dp, beta)
    case ('ls')
      beta = -dot_product(y, w) / &
             (dot_product(state%previous_direction, state%previous_gradient) + eps)
    case ('hz', 'hz+')
      if (canonical_name(state%control%preconditioner) == 'l-bfgs' .and. &
          state%history_count > 0) then
        call lbfgs_apply(state, y, v)
      else
        v = y
      end if
      denom = dot_product(state%previous_direction, y) + eps
      beta = dot_product(y, w) / denom - &
             2.0_dp * dot_product(y, v) * dot_product(state%previous_direction, &
             state%gradient) / (denom * denom)
      if (update == 'hz+') then
        eta_k = -1.0_dp / (max(dot_product(state%previous_direction, &
                state%previous_direction), eps) * &
                min(0.01_dp, max(vector_norm(state%previous_gradient), eps)))
        beta = max(eta_k, beta)
      end if
    case ('prfr')
      beta = dot_product(w, y) / &
             (dot_product(state%previous_gradient, wold) + eps)
      denom = dot_product(state%gradient, w) / &
              (dot_product(state%previous_gradient, wold) + eps)
      beta = max(-denom, min(denom, beta))
    case default
      beta = 0.0_dp
    end select

    direction = -w + beta * state%previous_direction
    if (dot_product(state%gradient, direction) >= 0.0_dp) direction = -w
  end subroutine cg_direction_value

  subroutine obtain_hessian(state, fg, h, hessian, user_data, ok)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(out) :: h(:, :)
    procedure(mize_hessian_interface), optional :: hessian
    class(*), intent(inout), optional :: user_data
    logical, intent(out) :: ok

    integer :: i
    real(dp), allocatable :: xp(:), gp(:)
    real(dp) :: fp, step
    logical :: eval_ok

    if (present(hessian)) then
      call invoke_hessian(hessian, state%par, h, user_data)
      ok = all(ieee_is_finite(h))
      return
    end if

    allocate(xp(state%n), gp(state%n))
    do i = 1, state%n
      xp = state%par
      step = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(state%par(i)))
      xp(i) = xp(i) + step
      call evaluate_fg(fg, xp, fp, gp, user_data, eval_ok)
      state%function_evaluations = state%function_evaluations + 1
      state%gradient_evaluations = state%gradient_evaluations + 1
      if (.not. eval_ok) then
        ok = .false.
        return
      end if
      h(:, i) = (gp - state%gradient) / step
    end do
    h = 0.5_dp * (h + transpose(h))
    ok = all(ieee_is_finite(h))
  end subroutine obtain_hessian

  subroutine truncated_newton_direction(state, fg, direction, hessian, hvp, user_data, ok)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(out) :: direction(:)
    procedure(mize_hessian_interface), optional :: hessian
    procedure(mize_hvp_interface), optional :: hvp
    class(*), intent(inout), optional :: user_data
    logical, intent(out) :: ok

    real(dp), allocatable :: z(:), r(:), p(:), hp(:), y(:), znew(:)
    real(dp) :: rr, rr_new, alpha, beta, pHp, tol, gn
    integer :: j
    logical :: hv_ok

    allocate(z(state%n), r(state%n), p(state%n), hp(state%n), y(state%n), znew(state%n))
    if (canonical_name(state%control%tn_init) == 'previous' .and. state%have_previous) then
      z = state%previous_direction
      call apply_hessian_vector(state, fg, z, hp, hessian, hvp, user_data, hv_ok)
      if (.not. hv_ok) then
        ok = .false.
        return
      end if
      r = state%gradient + hp
    else if (canonical_name(state%control%tn_init) == 'l-bfgs' .and. &
             state%history_count > 0) then
      call lbfgs_apply(state, state%gradient, z)
      z = -z
      call apply_hessian_vector(state, fg, z, hp, hessian, hvp, user_data, hv_ok)
      if (.not. hv_ok) then
        ok = .false.
        return
      end if
      r = state%gradient + hp
    else
      z = 0.0_dp
      r = state%gradient
    end if

    if (canonical_name(state%control%preconditioner) == 'l-bfgs' .and. &
        state%history_count > 0) then
      call lbfgs_apply(state, r, y)
    else
      y = r
    end if
    p = -y
    rr = dot_product(r, y)
    gn = vector_norm(state%gradient)
    tol = min(0.5_dp, sqrt(max(gn, 0.0_dp))) * gn
    ok = .true.

    if (.not. ieee_is_finite(rr) .or. rr <= 0.0_dp) then
      direction = -state%gradient
      return
    end if

    do j = 1, min(state%n, state%control%tn_max_iterations)
      call apply_hessian_vector(state, fg, p, hp, hessian, hvp, user_data, hv_ok)
      if (.not. hv_ok) then
        ok = .false.
        return
      end if
      pHp = dot_product(p, hp)
      if (.not. ieee_is_finite(pHp) .or. pHp <= 0.0_dp) then
        if (j == 1) z = -state%gradient
        exit
      end if
      alpha = rr / pHp
      znew = z + alpha * p
      r = r + alpha * hp
      if (vector_norm(r) <= tol) then
        z = znew
        exit
      end if
      if (canonical_name(state%control%preconditioner) == 'l-bfgs' .and. &
          state%history_count > 0) then
        call lbfgs_apply(state, r, y)
      else
        y = r
      end if
      rr_new = dot_product(r, y)
      if (.not. ieee_is_finite(rr_new) .or. rr_new <= 0.0_dp) exit
      beta = rr_new / rr
      p = -y + beta * p
      z = znew
      rr = rr_new
    end do
    direction = z
    if (dot_product(state%gradient, direction) >= 0.0_dp) direction = -state%gradient
  end subroutine truncated_newton_direction

  subroutine apply_hessian_vector(state, fg, v, hv, hessian, hvp, user_data, ok)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: v(:)
    real(dp), intent(out) :: hv(:)
    procedure(mize_hessian_interface), optional :: hessian
    procedure(mize_hvp_interface), optional :: hvp
    class(*), intent(inout), optional :: user_data
    logical, intent(out) :: ok

    real(dp), allocatable :: h(:, :), xp(:), gp(:)
    real(dp) :: fp, step, vn
    logical :: eval_ok

    if (present(hvp)) then
      call invoke_hvp(hvp, state%par, v, hv, user_data)
      ok = all(ieee_is_finite(hv))
      return
    end if
    if (present(hessian)) then
      allocate(h(state%n, state%n))
      call invoke_hessian(hessian, state%par, h, user_data)
      hv = matmul(h, v)
      ok = all(ieee_is_finite(hv))
      return
    end if

    vn = vector_norm(v)
    if (vn <= tiny(1.0_dp)) then
      hv = 0.0_dp
      ok = .true.
      return
    end if
    step = sqrt(epsilon(1.0_dp)) * (1.0_dp + vector_norm(state%par)) / vn
    allocate(xp(state%n), gp(state%n))
    xp = state%par + step * v
    call evaluate_fg(fg, xp, fp, gp, user_data, eval_ok)
    state%function_evaluations = state%function_evaluations + 1
    state%gradient_evaluations = state%gradient_evaluations + 1
    if (.not. eval_ok) then
      ok = .false.
      return
    end if
    hv = (gp - state%gradient) / step
    ok = all(ieee_is_finite(hv))
  end subroutine apply_hessian_vector

  subroutine dbd_step(state, fg, x_new, f_new, g_new, user_data, ok)
    type(mize_state_t), intent(inout) :: state
    procedure(mize_fg_interface) :: fg
    real(dp), intent(out) :: x_new(:), f_new, g_new(:)
    class(*), intent(inout), optional :: user_data
    logical, intent(out) :: ok

    integer :: i
    logical :: same_sign

    do i = 1, state%n
      if (state%iteration == 0 .or. abs(state%dbd_average(i)) <= tiny(1.0_dp)) then
        same_sign = .true.
      else
        same_sign = state%dbd_average(i) * state%gradient(i) >= 0.0_dp
      end if
      if (same_sign) then
        if (state%control%step_up_add) then
          state%dbd_scale(i) = state%dbd_scale(i) + state%control%step_up
        else
          state%dbd_scale(i) = state%dbd_scale(i) * state%control%step_up
        end if
      else
        state%dbd_scale(i) = state%dbd_scale(i) * state%control%step_down
      end if
      state%dbd_scale(i) = max(state%control%min_step, &
                               min(state%control%max_step, state%dbd_scale(i)))
    end do

    state%dbd_average = (1.0_dp - state%control%dbd_weight) * state%gradient + &
                        state%control%dbd_weight * state%dbd_average
    x_new = state%par - state%dbd_scale * state%gradient
    call evaluate_fg(fg, x_new, f_new, g_new, user_data, ok)
    state%function_evaluations = state%function_evaluations + 1
    state%gradient_evaluations = state%gradient_evaluations + 1
  end subroutine dbd_step

  real(dp) function momentum_value(state, momentum_schedule, user_data) result(mu)
    type(mize_state_t), intent(in) :: state
    procedure(mize_momentum_interface), optional :: momentum_schedule
    class(*), intent(inout), optional :: user_data
    character(len=:), allocatable :: schedule
    real(dp) :: t

    if (present(momentum_schedule)) then
      if (present(user_data)) then
        call momentum_schedule(state%iteration + 1, state%control%max_iterations, mu, user_data)
      else
        call momentum_schedule(state%iteration + 1, state%control%max_iterations, mu)
      end if
      mu = max(0.0_dp, min(1.0_dp, mu))
      return
    end if

    schedule = canonical_name(state%control%mom_schedule)
    select case (schedule)
    case ('ramp')
      t = real(state%iteration, dp) / real(max(1, state%control%max_iterations - 1), dp)
      mu = state%control%mom_init + t * (state%control%mom_final - state%control%mom_init)
    case ('switch')
      if (state%iteration < state%control%mom_switch_iter) then
        mu = state%control%mom_init
      else
        mu = state%control%mom_final
      end if
    case ('nsconvex', 'nesterov')
      mu = nesterov_momentum(state%iteration + 1, state%control%nest_q, &
                             state%control%nest_burn_in)
    case default
      mu = state%control%mom_init
    end select
    mu = max(0.0_dp, min(1.0_dp, mu))
  end function momentum_value

  real(dp) function nesterov_momentum(iteration, q, burn_in) result(mu)
    integer, intent(in) :: iteration, burn_in
    real(dp), intent(in) :: q
    real(dp) :: theta_old, theta_new, disc
    integer :: k

    if (iteration <= burn_in + 1) then
      mu = 0.0_dp
      return
    end if
    mu = 0.0_dp
    theta_old = 1.0_dp
    do k = 1, iteration - burn_in - 1
      disc = (theta_old * theta_old - q) ** 2 + &
             4.0_dp * theta_old * theta_old
      theta_new = 0.5_dp * (-(theta_old * theta_old - q) + sqrt(max(0.0_dp, disc)))
      if (theta_old <= tiny(1.0_dp)) exit
      mu = theta_old * (1.0_dp - theta_old) / &
           (theta_old * theta_old + theta_new)
      theta_old = theta_new
    end do
    mu = max(0.0_dp, min(1.0_dp, mu))
  end function nesterov_momentum

  subroutine apply_adaptive_restart(state, method, f_old, f_new, g_new, update)
    type(mize_state_t), intent(inout) :: state
    character(len=*), intent(in) :: method
    real(dp), intent(in) :: f_old, f_new, g_new(:), update(:)
    character(len=:), allocatable :: restart
    logical :: restart_now

    if (method /= 'momentum' .and. method /= 'nag') return
    restart = canonical_name(state%control%restart)
    if (len(restart) == 0) return
    if (state%iteration - state%last_restart <= state%control%restart_wait) return

    restart_now = .false.
    select case (restart)
    case ('fn', 'function')
      restart_now = f_new > f_old
    case ('gr', 'gradient')
      restart_now = dot_product(g_new, update) > 0.0_dp
    end select
    if (restart_now) then
      state%velocity = 0.0_dp
      state%last_restart = state%iteration
    end if
  end subroutine apply_adaptive_restart

  subroutine check_mize_convergence(state)
    type(mize_state_t), intent(inout) :: state
    real(dp) :: fn_change, fn_scale

    if (gradient_converged(state)) then
      call terminate_state(state, mize_success, 'gradient tolerance satisfied')
      return
    end if

    if (state%control%step_tol >= 0.0_dp) then
      if (state%last_step_norm <= state%control%step_tol * &
          (1.0_dp + vector_norm(state%par))) then
        call terminate_state(state, mize_success, 'step tolerance satisfied')
        return
      end if
    end if

    fn_change = abs(state%value - state%previous_value)
    fn_scale = max(abs(state%value), abs(state%previous_value))
    if (state%control%abs_tol >= 0.0_dp .and. fn_change <= state%control%abs_tol) then
      call terminate_state(state, mize_success, 'absolute function tolerance satisfied')
      return
    end if
    if (state%control%rel_tol >= 0.0_dp .and. &
        fn_change <= state%control%rel_tol * max(1.0_dp, fn_scale)) then
      call terminate_state(state, mize_success, 'relative function tolerance satisfied')
      return
    end if

    if (state%iteration >= state%control%max_iterations) then
      call terminate_state(state, mize_max_iterations, 'maximum iterations reached')
      return
    end if
    if (state%function_evaluations >= state%control%max_evaluations) then
      call terminate_state(state, mize_max_evaluations, 'maximum evaluations reached')
    end if
  end subroutine check_mize_convergence

  logical function gradient_converged(state) result(converged)
    type(mize_state_t), intent(in) :: state
    real(dp) :: g2, gi

    converged = .false.
    g2 = vector_norm(state%gradient)
    gi = maxval(abs(state%gradient))
    if (state%control%grad_tol >= 0.0_dp .and. g2 <= state%control%grad_tol) then
      converged = .true.
    end if
    if (state%control%ginf_tol >= 0.0_dp .and. gi <= state%control%ginf_tol) then
      converged = .true.
    end if
  end function gradient_converged

  subroutine mize_state_result(state, result)
    type(mize_state_t), intent(in) :: state
    type(mize_result_t), intent(out) :: result
    integer :: last

    if (allocated(state%par)) then
      allocate(result%par(size(state%par)), result%gradient(size(state%gradient)))
      allocate(result%best_par(size(state%best_par)))
      result%par = state%par
      result%gradient = state%gradient
      result%best_par = state%best_par
      result%value = state%value
      result%best_value = state%best_value
      result%gradient_norm = vector_norm(state%gradient)
      result%gradient_inf_norm = maxval(abs(state%gradient))
    end if
    result%iterations = state%iteration
    result%function_evaluations = state%function_evaluations
    result%gradient_evaluations = state%gradient_evaluations
    result%status = state%status
    result%converged = state%status == mize_success .and. state%terminated
    if (allocated(state%message)) then
      result%message = state%message
    else
      result%message = mize_status_message(state%status)
    end if

    if (state%control%store_progress .and. allocated(state%progress_value)) then
      last = min(state%iteration, ubound(state%progress_value, 1))
      allocate(result%progress_value(0:last))
      allocate(result%progress_gradient_norm(0:last))
      allocate(result%progress_step_norm(0:last))
      result%progress_value = state%progress_value(0:last)
      result%progress_gradient_norm = state%progress_gradient_norm(0:last)
      result%progress_step_norm = state%progress_step_norm(0:last)
    end if
  end subroutine mize_state_result

  subroutine check_mize_gradient(par, fg, check, method, relative_step, absolute_step, &
                                 tolerance, user_data)
    real(dp), intent(in) :: par(:)
    procedure(mize_fg_interface) :: fg
    type(gradient_check_t), intent(out) :: check
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: relative_step, absolute_step, tolerance
    class(*), intent(inout), optional :: user_data

    real(dp), allocatable :: xp(:), xm(:), gp(:), gm(:)
    real(dp) :: f0, fp, fm, rel, abs_step, tol, step, scale
    integer :: i, n
    logical :: ok
    character(len=:), allocatable :: meth

    n = size(par)
    allocate(check%analytic(n), check%finite_difference(n))
    allocate(check%absolute_error(n), check%relative_error(n))
    allocate(xp(n), xm(n), gp(n), gm(n))
    call evaluate_fg(fg, par, f0, check%analytic, user_data, ok)
    if (.not. ok) then
      check%passed = .false.
      return
    end if

    meth = 'central'
    if (present(method)) meth = canonical_name(method)
    rel = sqrt(epsilon(1.0_dp))
    if (present(relative_step)) rel = relative_step
    abs_step = 0.0_dp
    if (present(absolute_step)) abs_step = absolute_step
    tol = 1.0e-5_dp
    if (present(tolerance)) tol = tolerance

    do i = 1, n
      step = max(abs_step, rel * max(1.0_dp, abs(par(i))))
      xp = par
      xp(i) = xp(i) + step
      call evaluate_fg(fg, xp, fp, gp, user_data, ok)
      if (.not. ok) then
        check%finite_difference(i) = huge(1.0_dp)
        cycle
      end if
      if (meth == 'forward') then
        check%finite_difference(i) = (fp - f0) / step
      else
        xm = par
        xm(i) = xm(i) - step
        call evaluate_fg(fg, xm, fm, gm, user_data, ok)
        if (.not. ok) then
          check%finite_difference(i) = huge(1.0_dp)
          cycle
        end if
        check%finite_difference(i) = (fp - fm) / (2.0_dp * step)
      end if
    end do

    check%absolute_error = abs(check%analytic - check%finite_difference)
    do i = 1, n
      scale = max(1.0_dp, abs(check%analytic(i)), abs(check%finite_difference(i)))
      check%relative_error(i) = check%absolute_error(i) / scale
    end do
    check%maximum_absolute_error = maxval(check%absolute_error)
    check%maximum_relative_error = maxval(check%relative_error)
    check%worst_coordinate = maxloc(check%relative_error, dim=1)
    check%passed = check%maximum_relative_error <= tol
  end subroutine check_mize_gradient

  subroutine solve_linear_system(a_in, b_in, x, ok)
    real(dp), intent(in) :: a_in(:, :), b_in(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: a(:, :), b(:), row(:)
    real(dp) :: factor, pivot
    integer :: n, i, j, k, p

    n = size(b_in)
    allocate(a(n, n), b(n), row(n))
    a = a_in
    b = b_in
    ok = .true.

    do k = 1, n - 1
      p = k - 1 + maxloc(abs(a(k:n, k)), dim=1)
      pivot = a(p, k)
      if (.not. ieee_is_finite(pivot) .or. abs(pivot) <= 100.0_dp * epsilon(1.0_dp)) then
        ok = .false.
        x = 0.0_dp
        return
      end if
      if (p /= k) then
        row = a(k, :)
        a(k, :) = a(p, :)
        a(p, :) = row
        factor = b(k)
        b(k) = b(p)
        b(p) = factor
      end if
      do i = k + 1, n
        factor = a(i, k) / a(k, k)
        a(i, k:n) = a(i, k:n) - factor * a(k, k:n)
        b(i) = b(i) - factor * b(k)
      end do
    end do
    if (abs(a(n, n)) <= 100.0_dp * epsilon(1.0_dp)) then
      ok = .false.
      x = 0.0_dp
      return
    end if

    x = 0.0_dp
    do i = n, 1, -1
      x(i) = b(i)
      do j = i + 1, n
        x(i) = x(i) - a(i, j) * x(j)
      end do
      if (abs(a(i, i)) <= 100.0_dp * epsilon(1.0_dp)) then
        ok = .false.
        x = 0.0_dp
        return
      end if
      x(i) = x(i) / a(i, i)
    end do
    ok = all(ieee_is_finite(x))
  end subroutine solve_linear_system

  subroutine evaluate_fg(fg, x, f, g, user_data, ok)
    procedure(mize_fg_interface) :: fg
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data
    logical, intent(out) :: ok

    if (present(user_data)) then
      call fg(x, f, g, user_data)
    else
      call fg(x, f, g)
    end if
    ok = ieee_is_finite(f) .and. all(ieee_is_finite(g))
  end subroutine evaluate_fg

  subroutine invoke_hessian(hessian, x, h, user_data)
    procedure(mize_hessian_interface) :: hessian
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    class(*), intent(inout), optional :: user_data

    if (present(user_data)) then
      call hessian(x, h, user_data)
    else
      call hessian(x, h)
    end if
  end subroutine invoke_hessian

  subroutine invoke_hvp(hvp, x, v, hv, user_data)
    procedure(mize_hvp_interface) :: hvp
    real(dp), intent(in) :: x(:), v(:)
    real(dp), intent(out) :: hv(:)
    class(*), intent(inout), optional :: user_data

    if (present(user_data)) then
      call hvp(x, v, hv, user_data)
    else
      call hvp(x, v, hv)
    end if
  end subroutine invoke_hvp

  subroutine invoke_monitor(monitor, x, f, g, iteration, function_evaluations, &
                            gradient_evaluations, stop, user_data)
    procedure(mize_monitor_interface) :: monitor
    real(dp), intent(in) :: x(:), f, g(:)
    integer, intent(in) :: iteration, function_evaluations, gradient_evaluations
    logical, intent(out) :: stop
    class(*), intent(inout), optional :: user_data

    if (present(user_data)) then
      call monitor(x, f, g, iteration, function_evaluations, gradient_evaluations, &
                   stop, user_data)
    else
      call monitor(x, f, g, iteration, function_evaluations, gradient_evaluations, stop)
    end if
  end subroutine invoke_monitor

  pure real(dp) function vector_norm(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vector_norm

  subroutine normalize_vector(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: value
    value = vector_norm(x)
    if (value > 0.0_dp) x = x / value
  end subroutine normalize_vector

  pure logical function curvature_ok_update(y, s) result(ok)
    real(dp), intent(in) :: y(:), s(:)
    real(dp) :: ys, scale
    ys = dot_product(y, s)
    scale = 1.0e-8_dp * vector_norm(y) * vector_norm(s)
    ok = ieee_is_finite(ys) .and. ieee_is_finite(scale) .and. ys > scale
  end function curvature_ok_update

  real(dp) function initial_step(state) result(alpha)
    type(mize_state_t), intent(in) :: state
    character(len=:), allocatable :: method
    method = canonical_name(state%control%method)
    alpha = state%control%step0
    if (state%iteration > 0 .and. state%last_step > 0.0_dp) then
      if (method == 'bfgs' .or. method == 'l-bfgs' .or. method == 'newton' .or. &
          method == 'phess' .or. method == 'tn') then
        alpha = min(1.0_dp, max(state%control%min_step, state%last_step * 2.0_dp))
      else
        alpha = min(state%control%max_step, max(state%control%min_step, &
                    state%last_step * 1.5_dp))
      end if
    end if
    alpha = max(state%control%min_step, min(state%control%max_step, alpha))
  end function initial_step

  logical function valid_control(control) result(ok)
    type(mize_control_t), intent(in) :: control
    character(len=:), allocatable :: method

    method = canonical_name(control%method)
    ok = any(method == [character(len=10) :: 'sd', 'bfgs', 'sr1', 'l-bfgs', 'cg', &
              'newton', 'phess', 'tn', 'nag', 'momentum', 'dbd'])
    ok = ok .and. control%memory >= 1
    ok = ok .and. control%max_iterations >= 0
    ok = ok .and. control%max_evaluations >= 1
    ok = ok .and. control%max_line_search >= 1
    ok = ok .and. control%step0 > 0.0_dp
    ok = ok .and. control%c1 > 0.0_dp .and. control%c1 < 1.0_dp
    ok = ok .and. control%c2 > control%c1 .and. control%c2 < 1.0_dp
    ok = ok .and. control%step_down > 0.0_dp .and. control%step_down < 1.0_dp
    ok = ok .and. control%step_up > 0.0_dp
    ok = ok .and. control%dbd_weight >= 0.0_dp .and. control%dbd_weight <= 1.0_dp
    ok = ok .and. control%nest_q >= 0.0_dp .and. control%nest_q <= 1.0_dp
  end function valid_control

  pure function canonical_name(name) result(lower)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: lower
    integer :: i, code

    lower = trim(adjustl(name))
    do i = 1, len(lower)
      code = iachar(lower(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lower(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
      if (lower(i:i) == '_') lower(i:i) = '-'
    end do
    select case (lower)
    case ('lbfgs')
      lower = 'l-bfgs'
    case ('mom')
      lower = 'momentum'
    case ('mt')
      lower = 'more-thuente'
    case ('hz')
      lower = 'hager-zhang'
    case ('minfunc')
      lower = 'schmidt'
    end select
  end function canonical_name

  subroutine fail_state(state, status, message)
    type(mize_state_t), intent(inout) :: state
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    state%status = status
    state%message = message
    state%initialized = .false.
    state%terminated = .true.
  end subroutine fail_state

  subroutine terminate_state(state, status, message)
    type(mize_state_t), intent(inout) :: state
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    state%status = status
    state%message = message
    state%terminated = .true.
  end subroutine terminate_state

  function mize_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message
    select case (status)
    case (mize_success)
      message = 'success'
    case (mize_max_iterations)
      message = 'maximum iterations reached'
    case (mize_max_evaluations)
      message = 'maximum evaluations reached'
    case (mize_line_search_failure)
      message = 'line search failed'
    case (mize_invalid_input)
      message = 'invalid input'
    case (mize_nonfinite)
      message = 'non-finite objective, gradient, or search direction'
    case (mize_user_stop)
      message = 'stopped by user callback'
    case default
      message = 'unknown status'
    end select
  end function mize_status_message

end module mize_mod
