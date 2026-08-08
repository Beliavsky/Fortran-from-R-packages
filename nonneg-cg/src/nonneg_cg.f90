module nonneg_cg
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    integer, parameter, public :: nonneg_cg_tol_achieved = 0
    integer, parameter, public :: nonneg_cg_stop_maxnfeval = 1
    integer, parameter, public :: nonneg_cg_stop_maxiter = 2
    integer, parameter, public :: nonneg_cg_out_of_memory = 3
    integer, parameter, public :: nonneg_cg_invalid_input = -1
    integer, parameter, public :: nonneg_cg_nonfinite_value = -2
    integer, parameter, public :: nonneg_cg_cancelled = -3

    type, public :: nonneg_cg_control_t
        real(dp) :: tol = 1.0e-4_dp
        integer :: maxnfeval = 1500
        integer :: maxiter = 200
        real(dp) :: decr_lnsrch = 0.5_dp
        real(dp) :: lnsrch_const = 0.01_dp
        integer :: max_ls = 20
        logical :: extra_nonneg_tol = .false.
        logical :: verbose = .false.
    end type nonneg_cg_control_t

    type, public :: nonneg_cg_result_t
        real(dp), allocatable :: x(:)
        real(dp) :: fun = huge(1.0_dp)
        integer :: niter = 0
        integer :: nfeval = 0
        integer :: objective_calls = 0
        integer :: status = nonneg_cg_invalid_input
        real(dp) :: directional_derivative = huge(1.0_dp)
        character(len=96) :: message = ''
    end type nonneg_cg_result_t

    abstract interface
        function objective_callback(x, user_data) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            class(*), intent(inout), optional :: user_data
            real(dp) :: value
        end function objective_callback

        subroutine gradient_callback(x, gradient, user_data)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: gradient(:)
            class(*), intent(inout), optional :: user_data
        end subroutine gradient_callback

        function monitor_callback(x, value, iteration, user_data) result(cancel)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: value
            integer, intent(in) :: iteration
            class(*), intent(inout), optional :: user_data
            logical :: cancel
        end function monitor_callback
    end interface

    public :: objective_callback
    public :: gradient_callback
    public :: monitor_callback
    public :: minimize_nonneg_cg
    public :: nonneg_cg_status_message

contains

    subroutine minimize_nonneg_cg(x, objective, gradient, result, control, user_data, monitor)
        real(dp), intent(inout) :: x(:)
        procedure(objective_callback) :: objective
        procedure(gradient_callback) :: gradient
        type(nonneg_cg_result_t), intent(out) :: result
        type(nonneg_cg_control_t), intent(in), optional :: control
        class(*), intent(inout), optional :: user_data
        procedure(monitor_callback), optional :: monitor

        type(nonneg_cg_control_t) :: ctrl
        real(dp), allocatable, target :: grad_store(:, :)
        real(dp), allocatable, target :: direction_store(:, :)
        real(dp), pointer :: grad_curr(:)
        real(dp), pointer :: grad_prev(:)
        real(dp), pointer :: direction_curr(:)
        real(dp), pointer :: direction_prev(:)
        real(dp) :: beta
        real(dp) :: curr_fun_val
        real(dp) :: direction_norm_sq
        real(dp) :: grad_prev_norm_sq
        real(dp) :: max_step
        real(dp) :: new_fun_val
        real(dp) :: prod_grad_dir
        real(dp) :: step_change
        real(dp) :: theta
        integer :: curr_slot
        integer :: i
        integer :: iter
        integer :: ls
        integer :: max_eval_internal
        integer :: max_iter_internal
        integer :: stat
        logical :: cancel
        logical :: revert_x

        ctrl = nonneg_cg_control_t()
        if (present(control)) ctrl = control

        call initialize_result(result, x)
        if (.not. valid_inputs(x, ctrl, result)) return

        max_eval_internal = ctrl%maxnfeval
        if (max_eval_internal == 0) max_eval_internal = huge(0)
        max_iter_internal = ctrl%maxiter
        if (max_iter_internal == 0) max_iter_internal = huge(0)

        allocate(grad_store(size(x), 2), direction_store(size(x), 2), stat=stat)
        if (stat /= 0) then
            call finish_result(result, x, nonneg_cg_out_of_memory, huge(1.0_dp), 0, 0, 0)
            return
        end if

        curr_fun_val = call_objective(objective, x, user_data)
        result%objective_calls = 1
        result%nfeval = 1
        if (.not. ieee_is_finite(curr_fun_val)) then
            call finish_result(result, x, nonneg_cg_nonfinite_value, curr_fun_val, 0, 1, 1)
            return
        end if

        if (ctrl%verbose) call print_header(size(x), curr_fun_val)

        curr_slot = 1
        revert_x = .false.
        max_step = 0.0_dp
        ls = 0
        grad_prev_norm_sq = 0.0_dp
        prod_grad_dir = huge(1.0_dp)

        do iter = 0, max_iter_internal - 1
            grad_curr => grad_store(:, curr_slot)
            direction_curr => direction_store(:, curr_slot)

            call call_gradient(gradient, x, grad_curr, user_data)
            if (any(.not. ieee_is_finite(grad_curr))) then
                call finish_result(result, x, nonneg_cg_nonfinite_value, curr_fun_val, iter, &
                    result%nfeval, result%objective_calls, prod_grad_dir)
                return
            end if

            do i = 1, size(x)
                if (x(i) <= 0.0_dp .and. grad_curr(i) >= 0.0_dp) then
                    direction_curr(i) = 0.0_dp
                else
                    direction_curr(i) = -grad_curr(i)
                end if
            end do

            if (iter > 0) then
                direction_prev => direction_store(:, 3 - curr_slot)
                grad_prev => grad_store(:, 3 - curr_slot)

                if (grad_prev_norm_sq <= tiny(1.0_dp)) then
                    prod_grad_dir = dot_vec(grad_curr, direction_curr)
                    call finish_result(result, x, nonneg_cg_tol_achieved, curr_fun_val, iter, &
                        result%nfeval, result%objective_calls, prod_grad_dir)
                    return
                end if

                theta = 0.0_dp
                beta = 0.0_dp
                do i = 1, size(x)
                    if (x(i) > 0.0_dp) then
                        theta = theta + grad_curr(i) * direction_prev(i)
                        beta = beta + grad_curr(i) * (grad_curr(i) - grad_prev(i))
                    end if
                end do
                theta = theta / grad_prev_norm_sq
                beta = beta / grad_prev_norm_sq

                do i = 1, size(x)
                    if (x(i) > 0.0_dp) then
                        direction_curr(i) = direction_curr(i) + beta * direction_prev(i) - &
                            theta * (grad_curr(i) - grad_prev(i))
                    end if
                end do
            end if

            prod_grad_dir = dot_vec(grad_curr, direction_curr)
            result%directional_derivative = prod_grad_dir
            if (abs(prod_grad_dir) <= ctrl%tol) then
                call finish_result(result, x, nonneg_cg_tol_achieved, curr_fun_val, iter, &
                    result%nfeval, result%objective_calls, prod_grad_dir)
                return
            end if

            max_step = 1.0_dp
            do i = 1, size(x)
                if (direction_curr(i) < 0.0_dp) then
                    max_step = min(max_step, -x(i) / direction_curr(i))
                end if
            end do

            call axpy_inplace(max_step, direction_curr, x)
            direction_norm_sq = dot_vec(direction_curr, direction_curr)

            do ls = 0, ctrl%max_ls - 1
                if (ctrl%extra_nonneg_tol) then
                    where (x <= 0.0_dp) x = 0.0_dp
                end if

                new_fun_val = call_objective(objective, x, user_data)
                result%objective_calls = result%objective_calls + 1

                if (ieee_is_finite(new_fun_val)) then
                    if (new_fun_val <= curr_fun_val - ctrl%lnsrch_const * &
                        (max_step * ctrl%decr_lnsrch**real(ls, dp))**2 * direction_norm_sq) exit
                end if

                result%nfeval = result%nfeval + 1
                if (result%nfeval >= max_eval_internal) then
                    revert_x = .true.
                    exit
                end if

                step_change = max_step * (ctrl%decr_lnsrch**real(ls + 1, dp) - ctrl%decr_lnsrch**real(ls, dp))
                call axpy_inplace(step_change, direction_curr, x)
            end do

            if (revert_x) then
                call axpy_inplace(-max_step * ctrl%decr_lnsrch**real(ls, dp), direction_curr, x)
                if (ctrl%extra_nonneg_tol) then
                    where (x <= 0.0_dp) x = 0.0_dp
                end if
                call finish_result(result, x, nonneg_cg_stop_maxnfeval, curr_fun_val, iter, &
                    result%nfeval, result%objective_calls, prod_grad_dir)
                if (ctrl%verbose) call print_footer(result)
                return
            end if

            curr_fun_val = new_fun_val
            if (present(monitor)) then
                cancel = call_monitor(monitor, x, curr_fun_val, iter + 1, user_data)
                if (cancel) then
                    call finish_result(result, x, nonneg_cg_cancelled, curr_fun_val, iter + 1, &
                        result%nfeval, result%objective_calls, prod_grad_dir)
                    if (ctrl%verbose) call print_footer(result)
                    return
                end if
            end if

            grad_prev_norm_sq = dot_vec(grad_curr, grad_curr)
            curr_slot = 3 - curr_slot

            if (ctrl%verbose) then
                write (*, '(a,i0,a,es14.6,a,es14.6,a,i0,a,i0)') &
                    'iteration ', iter + 1, ': f = ', curr_fun_val, ', |g.d| = ', &
                    abs(prod_grad_dir), ', nfeval = ', result%nfeval, ', ls = ', ls + 1
            end if
        end do

        call finish_result(result, x, nonneg_cg_stop_maxiter, curr_fun_val, max_iter_internal, &
            result%nfeval, result%objective_calls, prod_grad_dir)
        if (ctrl%verbose) call print_footer(result)
    end subroutine minimize_nonneg_cg

    function dot_vec(x, y) result(value)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: y(:)
        real(dp) :: value
        integer :: i

        value = 0.0_dp
        do i = 1, size(x)
            value = value + x(i) * y(i)
        end do
    end function dot_vec

    subroutine axpy_inplace(alpha, x, y)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: x(:)
        real(dp), intent(inout) :: y(:)
        integer :: i

        do i = 1, size(x)
            y(i) = y(i) + alpha * x(i)
        end do
    end subroutine axpy_inplace

    function call_objective(objective, x, user_data) result(value)
        procedure(objective_callback) :: objective
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: value

        if (present(user_data)) then
            value = objective(x, user_data)
        else
            value = objective(x)
        end if
    end function call_objective

    subroutine call_gradient(gradient, x, grad, user_data)
        procedure(gradient_callback) :: gradient
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: grad(:)
        class(*), intent(inout), optional :: user_data

        if (present(user_data)) then
            call gradient(x, grad, user_data)
        else
            call gradient(x, grad)
        end if
    end subroutine call_gradient

    function call_monitor(monitor, x, value, iteration, user_data) result(cancel)
        procedure(monitor_callback) :: monitor
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: value
        integer, intent(in) :: iteration
        class(*), intent(inout), optional :: user_data
        logical :: cancel

        if (present(user_data)) then
            cancel = monitor(x, value, iteration, user_data)
        else
            cancel = monitor(x, value, iteration)
        end if
    end function call_monitor

    function valid_inputs(x, control, result) result(ok)
        real(dp), intent(in) :: x(:)
        type(nonneg_cg_control_t), intent(in) :: control
        type(nonneg_cg_result_t), intent(inout) :: result
        logical :: ok

        ok = .false.
        if (size(x) == 0) then
            result%message = 'x must contain at least one variable'
        else if (any(.not. ieee_is_finite(x))) then
            result%message = 'x must contain only finite values'
        else if (any(x < 0.0_dp)) then
            result%message = 'x must be feasible: every component must be nonnegative'
        else if (control%tol < 0.0_dp .or. .not. ieee_is_finite(control%tol)) then
            result%message = 'tol must be finite and nonnegative'
        else if (control%maxnfeval < 0) then
            result%message = 'maxnfeval must be nonnegative; zero means unlimited'
        else if (control%maxiter < 0) then
            result%message = 'maxiter must be nonnegative; zero means unlimited'
        else if (control%decr_lnsrch <= 0.0_dp .or. control%decr_lnsrch >= 1.0_dp) then
            result%message = 'decr_lnsrch must lie strictly between zero and one'
        else if (control%lnsrch_const <= 0.0_dp .or. control%lnsrch_const >= 1.0_dp) then
            result%message = 'lnsrch_const must lie strictly between zero and one'
        else if (control%max_ls <= 0) then
            result%message = 'max_ls must be positive'
        else
            ok = .true.
            return
        end if
        result%status = nonneg_cg_invalid_input
    end function valid_inputs

    subroutine initialize_result(result, x)
        type(nonneg_cg_result_t), intent(out) :: result
        real(dp), intent(in) :: x(:)

        allocate(result%x(size(x)))
        result%x = x
        result%fun = huge(1.0_dp)
        result%niter = 0
        result%nfeval = 0
        result%objective_calls = 0
        result%status = nonneg_cg_invalid_input
        result%directional_derivative = huge(1.0_dp)
        result%message = ''
    end subroutine initialize_result

    subroutine finish_result(result, x, status, value, niter, nfeval, objective_calls, directional_derivative)
        type(nonneg_cg_result_t), intent(inout) :: result
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: status
        real(dp), intent(in) :: value
        integer, intent(in) :: niter
        integer, intent(in) :: nfeval
        integer, intent(in) :: objective_calls
        real(dp), intent(in), optional :: directional_derivative

        result%x = x
        result%fun = value
        result%niter = niter
        result%nfeval = nfeval
        result%objective_calls = objective_calls
        result%status = status
        if (present(directional_derivative)) result%directional_derivative = directional_derivative
        result%message = nonneg_cg_status_message(status)
    end subroutine finish_result

    function nonneg_cg_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=96) :: message

        select case (status)
        case (nonneg_cg_tol_achieved)
            message = 'tolerance achieved'
        case (nonneg_cg_stop_maxnfeval)
            message = 'maximum reported function evaluations reached'
        case (nonneg_cg_stop_maxiter)
            message = 'maximum conjugate-gradient iterations reached'
        case (nonneg_cg_out_of_memory)
            message = 'memory allocation failed'
        case (nonneg_cg_invalid_input)
            message = 'invalid input'
        case (nonneg_cg_nonfinite_value)
            message = 'objective or gradient returned a non-finite value'
        case (nonneg_cg_cancelled)
            message = 'cancelled by monitor callback'
        case default
            message = 'unknown status'
        end select
    end function nonneg_cg_status_message

    subroutine print_header(n, value)
        integer, intent(in) :: n
        real(dp), intent(in) :: value

        write (*, '(a)') '********************************************'
        write (*, '(a)') 'Non-negative Conjugate Gradient Optimization'
        write (*, '(a,i0)') 'Number of variables: ', n
        write (*, '(a,es14.6)') 'Initial function value: ', value
    end subroutine print_header

    subroutine print_footer(result)
        type(nonneg_cg_result_t), intent(in) :: result

        write (*, '(a,a)') 'Terminated: ', trim(result%message)
        write (*, '(a,es14.6)') 'Last f(x): ', result%fun
    end subroutine print_footer

end module nonneg_cg
