module lbfgs_solver
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use lbfgs_kinds, only : dp
    use lbfgs_status
    implicit none
    private

    type, public :: lbfgs_parameter_t
        integer :: m = 6
        real(dp) :: epsilon = 1.0e-5_dp
        integer :: past = 0
        real(dp) :: delta = 0.0_dp
        integer :: max_iterations = 0
        integer :: linesearch = lbfgs_linesearch_default
        integer :: max_linesearch = 20
        real(dp) :: min_step = 1.0e-20_dp
        real(dp) :: max_step = 1.0e20_dp
        real(dp) :: ftol = 1.0e-4_dp
        real(dp) :: wolfe = 0.9_dp
        real(dp) :: gtol = 0.9_dp
        real(dp) :: xtol = epsilon(1.0_dp)
        real(dp) :: orthantwise_c = 0.0_dp
        integer :: orthantwise_start = 1
        integer :: orthantwise_end = 0
    contains
        procedure :: reset => reset_parameters
    end type lbfgs_parameter_t

    type, public :: lbfgs_result_t
        real(dp) :: value = huge(1.0_dp)
        integer :: status = lbfgserr_unknownerror
        integer :: iterations = 0
        integer :: evaluations = 0
        integer :: last_linesearch_evaluations = 0
        real(dp) :: x_norm = 0.0_dp
        real(dp) :: gradient_norm = 0.0_dp
        character(len=:), allocatable :: message
    end type lbfgs_result_t

    abstract interface
        subroutine lbfgs_evaluate_proc(x, f, g, step, user_data)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f
            real(dp), intent(out) :: g(:)
            real(dp), intent(in) :: step
            class(*), intent(inout), optional :: user_data
        end subroutine lbfgs_evaluate_proc

        subroutine lbfgs_objective_proc(x, f, user_data)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f
            class(*), intent(inout), optional :: user_data
        end subroutine lbfgs_objective_proc

        subroutine lbfgs_gradient_proc(x, g, user_data)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: g(:)
            class(*), intent(inout), optional :: user_data
        end subroutine lbfgs_gradient_proc

        integer function lbfgs_progress_proc(x, g, f, xnorm, gnorm, step, &
                iteration, line_evaluations, user_data)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: g(:)
            real(dp), intent(in) :: f
            real(dp), intent(in) :: xnorm
            real(dp), intent(in) :: gnorm
            real(dp), intent(in) :: step
            integer, intent(in) :: iteration
            integer, intent(in) :: line_evaluations
            class(*), intent(inout), optional :: user_data
        end function lbfgs_progress_proc
    end interface

    interface lbfgs_minimize
        module procedure lbfgs_minimize_combined
        module procedure lbfgs_minimize_separate
    end interface lbfgs_minimize

    public :: lbfgs_minimize
    public :: lbfgs_evaluate_proc, lbfgs_objective_proc, lbfgs_gradient_proc
    public :: lbfgs_progress_proc

    type :: iteration_data_t
        real(dp) :: alpha = 0.0_dp
        real(dp) :: ys = 0.0_dp
        real(dp), allocatable :: s(:)
        real(dp), allocatable :: y(:)
    end type iteration_data_t

contains

    subroutine reset_parameters(self)
        class(lbfgs_parameter_t), intent(inout) :: self

        self%m = 6
        self%epsilon = 1.0e-5_dp
        self%past = 0
        self%delta = 0.0_dp
        self%max_iterations = 0
        self%linesearch = lbfgs_linesearch_default
        self%max_linesearch = 20
        self%min_step = 1.0e-20_dp
        self%max_step = 1.0e20_dp
        self%ftol = 1.0e-4_dp
        self%wolfe = 0.9_dp
        self%gtol = 0.9_dp
        self%xtol = epsilon(1.0_dp)
        self%orthantwise_c = 0.0_dp
        self%orthantwise_start = 1
        self%orthantwise_end = 0
    end subroutine reset_parameters

    subroutine lbfgs_minimize_separate(objective, gradient, x, result, &
            parameters, progress, user_data)
        procedure(lbfgs_objective_proc) :: objective
        procedure(lbfgs_gradient_proc) :: gradient
        real(dp), intent(inout) :: x(:)
        type(lbfgs_result_t), intent(out) :: result
        type(lbfgs_parameter_t), intent(in), optional :: parameters
        procedure(lbfgs_progress_proc), optional :: progress
        class(*), intent(inout), optional :: user_data

        if (present(parameters)) then
            if (present(progress)) then
                if (present(user_data)) then
                    call lbfgs_minimize_combined(joined_evaluate, x, result, &
                        parameters, progress, user_data)
                else
                    call lbfgs_minimize_combined(joined_evaluate, x, result, &
                        parameters, progress)
                end if
            else
                if (present(user_data)) then
                    call lbfgs_minimize_combined(joined_evaluate, x, result, &
                        parameters, user_data=user_data)
                else
                    call lbfgs_minimize_combined(joined_evaluate, x, result, parameters)
                end if
            end if
        else
            if (present(progress)) then
                if (present(user_data)) then
                    call lbfgs_minimize_combined(joined_evaluate, x, result, &
                        progress=progress, user_data=user_data)
                else
                    call lbfgs_minimize_combined(joined_evaluate, x, result, progress=progress)
                end if
            else
                if (present(user_data)) then
                    call lbfgs_minimize_combined(joined_evaluate, x, result, &
                        user_data=user_data)
                else
                    call lbfgs_minimize_combined(joined_evaluate, x, result)
                end if
            end if
        end if

    contains

        subroutine joined_evaluate(x_current, f, g, step, ignored_data)
            real(dp), intent(in) :: x_current(:)
            real(dp), intent(out) :: f
            real(dp), intent(out) :: g(:)
            real(dp), intent(in) :: step
            class(*), intent(inout), optional :: ignored_data

            if (present(user_data)) then
                call objective(x_current, f, user_data)
                call gradient(x_current, g, user_data)
            else
                call objective(x_current, f)
                call gradient(x_current, g)
            end if
        end subroutine joined_evaluate

    end subroutine lbfgs_minimize_separate

    subroutine lbfgs_minimize_combined(evaluate, x, result, parameters, &
            progress, user_data)
        procedure(lbfgs_evaluate_proc) :: evaluate
        real(dp), intent(inout) :: x(:)
        type(lbfgs_result_t), intent(out) :: result
        type(lbfgs_parameter_t), intent(in), optional :: parameters
        procedure(lbfgs_progress_proc), optional :: progress
        class(*), intent(inout), optional :: user_data

        type(lbfgs_parameter_t) :: param
        type(iteration_data_t), allocatable :: history(:)
        real(dp), allocatable :: xp(:), g(:), gp(:), pg(:), d(:), work(:), pf(:)
        real(dp) :: f, step, xnorm, gnorm, ys, yy, beta, rate, denominator
        integer :: n, i, j, k, hist_end, bound, ls, status, start_index, end_index
        integer :: evaluations
        integer :: callback_status
        logical :: use_owlqn

        param = lbfgs_parameter_t()
        if (present(parameters)) param = parameters

        n = size(x)
        status = validate_parameters(n, param)
        if (status /= lbfgs_success) then
            call finish_result(result, status, huge(1.0_dp), 0, 0, 0, 0.0_dp, 0.0_dp)
            return
        end if

        start_index = param%orthantwise_start
        if (param%orthantwise_end == 0) then
            end_index = n
        else
            end_index = param%orthantwise_end
        end if
        use_owlqn = param%orthantwise_c /= 0.0_dp

        allocate(xp(n), g(n), gp(n), d(n), work(n))
        if (use_owlqn) allocate(pg(n))
        allocate(history(param%m))
        do i = 1, param%m
            allocate(history(i)%s(n), history(i)%y(n))
            history(i)%s = 0.0_dp
            history(i)%y = 0.0_dp
        end do
        if (param%past > 0) allocate(pf(param%past))

        evaluations = 0
        call evaluate_callback(evaluate, x, f, g, 0.0_dp, evaluations, user_data)
        if (.not. finite_evaluation(f, g)) then
            call finish_result(result, lbfgserr_logicerror, f, 0, evaluations, &
                0, vector_norm(x), vector_norm(g))
            return
        end if

        if (use_owlqn) then
            f = f + param%orthantwise_c * owlqn_x1norm(x, start_index, end_index)
            call owlqn_pseudo_gradient(pg, x, g, param%orthantwise_c, &
                start_index, end_index)
            d = -pg
        else
            d = -g
        end if

        if (allocated(pf)) pf(1) = f

        xnorm = vector_norm(x)
        if (use_owlqn) then
            gnorm = vector_norm(pg)
        else
            gnorm = vector_norm(g)
        end if
        if (gnorm / max(1.0_dp, xnorm) <= param%epsilon) then
            call finish_result(result, lbfgs_already_minimized, f, 0, &
                evaluations, 0, xnorm, gnorm)
            return
        end if

        step = 1.0_dp / max(vector_norm(d), tiny(1.0_dp))
        k = 1
        hist_end = 1
        ls = 0
        status = lbfgserr_unknownerror

        do
            xp = x
            gp = g

            if (use_owlqn) then
                ls = line_search_owlqn(evaluate, x, f, g, d, step, xp, pg, work, &
                    param, start_index, end_index, evaluations, user_data)
                if (ls >= 0) then
                    call owlqn_pseudo_gradient(pg, x, g, param%orthantwise_c, &
                        start_index, end_index)
                end if
            else if (param%linesearch == lbfgs_linesearch_morethuente) then
                ls = line_search_morethuente(evaluate, x, f, g, d, step, xp, &
                    param, evaluations, user_data)
            else
                ls = line_search_backtracking(evaluate, x, f, g, d, step, xp, &
                    param, evaluations, user_data)
            end if

            if (ls < 0) then
                x = xp
                g = gp
                status = ls
                exit
            end if

            xnorm = vector_norm(x)
            if (use_owlqn) then
                gnorm = vector_norm(pg)
            else
                gnorm = vector_norm(g)
            end if

            if (present(progress)) then
                if (present(user_data)) then
                    callback_status = progress(x, g, f, xnorm, gnorm, step, k, ls, user_data)
                else
                    callback_status = progress(x, g, f, xnorm, gnorm, step, k, ls)
                end if
                if (callback_status /= 0) then
                    status = lbfgserr_canceled
                    exit
                end if
            end if

            if (gnorm / max(1.0_dp, xnorm) <= param%epsilon) then
                status = lbfgs_success
                exit
            end if

            if (allocated(pf)) then
                if (k >= param%past) then
                    denominator = max(abs(f), 1.0_dp)
                    rate = (pf(mod(k, param%past) + 1) - f) / denominator
                    if (rate < param%delta) then
                        status = lbfgs_stop
                        exit
                    end if
                end if
                pf(mod(k, param%past) + 1) = f
            end if

            if (param%max_iterations > 0 .and. k >= param%max_iterations) then
                status = lbfgserr_maximumiteration
                exit
            end if

            history(hist_end)%s = x - xp
            history(hist_end)%y = g - gp
            ys = dot_product(history(hist_end)%y, history(hist_end)%s)
            yy = dot_product(history(hist_end)%y, history(hist_end)%y)
            history(hist_end)%ys = ys

            if (.not. ieee_is_finite(ys) .or. .not. ieee_is_finite(yy) .or. &
                    ys <= tiny(1.0_dp) .or. yy <= tiny(1.0_dp)) then
                if (use_owlqn) then
                    d = -pg
                else
                    d = -g
                end if
                step = 1.0_dp
                k = k + 1
                cycle
            end if

            bound = min(param%m, k)
            k = k + 1
            hist_end = mod(hist_end, param%m) + 1

            if (use_owlqn) then
                d = -pg
            else
                d = -g
            end if

            j = hist_end
            do i = 1, bound
                j = modulo(j - 2, param%m) + 1
                history(j)%alpha = dot_product(history(j)%s, d) / history(j)%ys
                d = d - history(j)%alpha * history(j)%y
            end do

            d = d * (ys / yy)

            do i = 1, bound
                beta = dot_product(history(j)%y, d) / history(j)%ys
                d = d + (history(j)%alpha - beta) * history(j)%s
                j = mod(j, param%m) + 1
            end do

            if (use_owlqn) then
                do i = start_index, end_index
                    if (d(i) * pg(i) >= 0.0_dp) d(i) = 0.0_dp
                end do
            end if

            step = 1.0_dp
        end do

        call finish_result(result, status, f, k, evaluations, ls, xnorm, gnorm)
    end subroutine lbfgs_minimize_combined

    integer function validate_parameters(n, param) result(status)
        integer, intent(in) :: n
        type(lbfgs_parameter_t), intent(in) :: param
        integer :: end_index

        status = lbfgs_success
        if (n <= 0 .or. param%m <= 0) then
            status = lbfgserr_invalid_n
        else if (param%epsilon < 0.0_dp) then
            status = lbfgserr_invalid_epsilon
        else if (param%past < 0) then
            status = lbfgserr_invalid_testperiod
        else if (param%delta < 0.0_dp) then
            status = lbfgserr_invalid_delta
        else if (param%min_step < 0.0_dp) then
            status = lbfgserr_invalid_minstep
        else if (param%max_step < param%min_step) then
            status = lbfgserr_invalid_maxstep
        else if (param%ftol < 0.0_dp .or. param%ftol >= 0.5_dp) then
            status = lbfgserr_invalid_ftol
        else if ((param%linesearch == lbfgs_linesearch_backtracking_wolfe .or. &
                param%linesearch == lbfgs_linesearch_backtracking_strong_wolfe) .and. &
                (param%wolfe <= param%ftol .or. param%wolfe >= 1.0_dp)) then
            status = lbfgserr_invalid_wolfe
        else if (param%gtol < 0.0_dp) then
            status = lbfgserr_invalid_gtol
        else if (param%xtol < 0.0_dp) then
            status = lbfgserr_invalid_xtol
        else if (param%max_linesearch <= 0) then
            status = lbfgserr_invalid_maxlinesearch
        else if (param%orthantwise_c < 0.0_dp) then
            status = lbfgserr_invalid_orthantwise
        else if (param%orthantwise_start < 1 .or. param%orthantwise_start > n) then
            status = lbfgserr_invalid_orthantwise_start
        else
            end_index = merge(n, param%orthantwise_end, param%orthantwise_end == 0)
            if (end_index < param%orthantwise_start .or. end_index > n) then
                status = lbfgserr_invalid_orthantwise_end
            else if (param%orthantwise_c /= 0.0_dp .and. &
                    param%linesearch /= lbfgs_linesearch_backtracking) then
                status = lbfgserr_invalid_linesearch
            else if (param%linesearch /= lbfgs_linesearch_morethuente .and. &
                    param%linesearch /= lbfgs_linesearch_backtracking_armijo .and. &
                    param%linesearch /= lbfgs_linesearch_backtracking_wolfe .and. &
                    param%linesearch /= lbfgs_linesearch_backtracking_strong_wolfe) then
                status = lbfgserr_invalid_linesearch
            end if
        end if
    end function validate_parameters

    subroutine evaluate_callback(evaluate, x, f, g, step, evaluations, user_data)
        procedure(lbfgs_evaluate_proc) :: evaluate
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(:)
        real(dp), intent(in) :: step
        integer, intent(inout) :: evaluations
        class(*), intent(inout), optional :: user_data

        if (present(user_data)) then
            call evaluate(x, f, g, step, user_data)
        else
            call evaluate(x, f, g, step)
        end if
        evaluations = evaluations + 1
    end subroutine evaluate_callback

    pure logical function finite_evaluation(f, g)
        real(dp), intent(in) :: f
        real(dp), intent(in) :: g(:)

        finite_evaluation = ieee_is_finite(f) .and. all(ieee_is_finite(g))
    end function finite_evaluation

    subroutine finish_result(result, status, value, iterations, evaluations, &
            line_evaluations, xnorm, gnorm)
        type(lbfgs_result_t), intent(out) :: result
        integer, intent(in) :: status, iterations, evaluations, line_evaluations
        real(dp), intent(in) :: value, xnorm, gnorm

        result%status = status
        result%value = value
        result%iterations = max(iterations, 0)
        result%evaluations = evaluations
        result%last_linesearch_evaluations = max(line_evaluations, 0)
        result%x_norm = xnorm
        result%gradient_norm = gnorm
        result%message = lbfgs_status_message(status)
    end subroutine finish_result

    pure real(dp) function vector_norm(x)
        real(dp), intent(in) :: x(:)

        vector_norm = sqrt(max(0.0_dp, dot_product(x, x)))
    end function vector_norm

    integer function line_search_backtracking(evaluate, x, f, g, direction, &
            step, xp, param, evaluations, user_data) result(count)
        procedure(lbfgs_evaluate_proc) :: evaluate
        real(dp), intent(inout) :: x(:), f, g(:), step
        real(dp), intent(in) :: direction(:), xp(:)
        type(lbfgs_parameter_t), intent(in) :: param
        integer, intent(inout) :: evaluations
        class(*), intent(inout), optional :: user_data

        real(dp), parameter :: decrease = 0.5_dp, increase = 2.1_dp
        real(dp) :: width, directional_gradient, initial_directional_gradient
        real(dp) :: initial_f, test_gradient

        count = 0
        if (step <= 0.0_dp) then
            count = lbfgserr_invalidparameters
            return
        end if

        initial_directional_gradient = dot_product(g, direction)
        if (initial_directional_gradient > 0.0_dp) then
            count = lbfgserr_increasegradient
            return
        end if

        initial_f = f
        test_gradient = param%ftol * initial_directional_gradient

        do
            x = xp + step * direction
            call evaluate_callback(evaluate, x, f, g, step, evaluations, user_data)
            count = count + 1

            if (.not. finite_evaluation(f, g)) then
                width = decrease
            else if (f > initial_f + step * test_gradient) then
                width = decrease
            else if (param%linesearch == lbfgs_linesearch_backtracking_armijo) then
                return
            else
                directional_gradient = dot_product(g, direction)
                if (directional_gradient < param%wolfe * initial_directional_gradient) then
                    width = increase
                else if (param%linesearch == lbfgs_linesearch_backtracking_wolfe) then
                    return
                else if (directional_gradient > -param%wolfe * &
                        initial_directional_gradient) then
                    width = decrease
                else
                    return
                end if
            end if

            if (step < param%min_step) then
                count = lbfgserr_minimumstep
                return
            else if (step > param%max_step) then
                count = lbfgserr_maximumstep
                return
            else if (count >= param%max_linesearch) then
                count = lbfgserr_maximumlinesearch
                return
            end if
            step = step * width
        end do
    end function line_search_backtracking

    integer function line_search_owlqn(evaluate, x, f, g, direction, step, xp, &
            pseudo_gradient, work, param, start_index, end_index, evaluations, &
            user_data) result(count)
        procedure(lbfgs_evaluate_proc) :: evaluate
        real(dp), intent(inout) :: x(:), f, g(:), step
        real(dp), intent(in) :: direction(:), xp(:), pseudo_gradient(:)
        real(dp), intent(out) :: work(:)
        type(lbfgs_parameter_t), intent(in) :: param
        integer, intent(in) :: start_index, end_index
        integer, intent(inout) :: evaluations
        class(*), intent(inout), optional :: user_data

        real(dp), parameter :: width = 0.5_dp
        real(dp) :: initial_f, test_gradient
        integer :: i

        count = 0
        if (step <= 0.0_dp) then
            count = lbfgserr_invalidparameters
            return
        end if

        do i = 1, size(x)
            if (xp(i) == 0.0_dp) then
                work(i) = -pseudo_gradient(i)
            else
                work(i) = xp(i)
            end if
        end do

        initial_f = f
        do
            x = xp + step * direction
            call owlqn_project(x, work, start_index, end_index)
            call evaluate_callback(evaluate, x, f, g, step, evaluations, user_data)
            if (finite_evaluation(f, g)) then
                f = f + param%orthantwise_c * owlqn_x1norm(x, start_index, end_index)
            else
                f = huge(1.0_dp)
            end if
            count = count + 1

            test_gradient = dot_product(x - xp, pseudo_gradient)
            if (f <= initial_f + param%ftol * test_gradient) return

            if (step < param%min_step) then
                count = lbfgserr_minimumstep
                return
            else if (step > param%max_step) then
                count = lbfgserr_maximumstep
                return
            else if (count >= param%max_linesearch) then
                count = lbfgserr_maximumlinesearch
                return
            end if
            step = step * width
        end do
    end function line_search_owlqn

    integer function line_search_morethuente(evaluate, x, f, g, direction, &
            step, xp, param, evaluations, user_data) result(count)
        procedure(lbfgs_evaluate_proc) :: evaluate
        real(dp), intent(inout) :: x(:), f, g(:), step
        real(dp), intent(in) :: direction(:), xp(:)
        type(lbfgs_parameter_t), intent(in) :: param
        integer, intent(inout) :: evaluations
        class(*), intent(inout), optional :: user_data

        logical :: bracketed, stage1
        integer :: update_status
        real(dp) :: dg, stx, fx, dgx, sty, fy, dgy
        real(dp) :: fxm, dgxm, fym, dgym, fm, dgm
        real(dp) :: initial_f, ftest1, dginit, dgtest
        real(dp) :: width, previous_width, stmin, stmax

        count = 0
        if (step <= 0.0_dp) then
            count = lbfgserr_invalidparameters
            return
        end if

        dginit = dot_product(g, direction)
        if (dginit > 0.0_dp) then
            count = lbfgserr_increasegradient
            return
        end if

        bracketed = .false.
        stage1 = .true.
        update_status = 0
        initial_f = f
        dgtest = param%ftol * dginit
        width = param%max_step - param%min_step
        previous_width = 2.0_dp * width
        stx = 0.0_dp
        sty = 0.0_dp
        fx = initial_f
        fy = initial_f
        dgx = dginit
        dgy = dginit

        do
            if (bracketed) then
                stmin = min(stx, sty)
                stmax = max(stx, sty)
            else
                stmin = stx
                stmax = step + 4.0_dp * (step - stx)
            end if

            step = max(param%min_step, min(param%max_step, step))

            if (bracketed .and. ((step <= stmin .or. step >= stmax) .or. &
                    count + 1 >= param%max_linesearch .or. update_status /= 0 .or. &
                    stmax - stmin <= param%xtol * stmax)) then
                step = stx
            end if

            x = xp + step * direction
            call evaluate_callback(evaluate, x, f, g, step, evaluations, user_data)
            count = count + 1
            if (.not. finite_evaluation(f, g)) then
                count = lbfgserr_rounding_error
                return
            end if
            dg = dot_product(g, direction)
            ftest1 = initial_f + step * dgtest

            if (bracketed .and. ((step <= stmin .or. step >= stmax) .or. &
                    update_status /= 0)) then
                count = lbfgserr_rounding_error
                return
            else if (step == param%max_step .and. f <= ftest1 .and. dg <= dgtest) then
                count = lbfgserr_maximumstep
                return
            else if (step == param%min_step .and. (f > ftest1 .or. dg >= dgtest)) then
                count = lbfgserr_minimumstep
                return
            else if (bracketed .and. stmax - stmin <= param%xtol * stmax) then
                count = lbfgserr_widthtoosmall
                return
            else if (count >= param%max_linesearch) then
                count = lbfgserr_maximumlinesearch
                return
            else if (f <= ftest1 .and. abs(dg) <= param%gtol * (-dginit)) then
                return
            end if

            if (stage1 .and. f <= ftest1 .and. &
                    min(param%ftol, param%gtol) * dginit <= dg) then
                stage1 = .false.
            end if

            if (stage1 .and. ftest1 < f .and. f <= fx) then
                fm = f - step * dgtest
                fxm = fx - stx * dgtest
                fym = fy - sty * dgtest
                dgm = dg - dgtest
                dgxm = dgx - dgtest
                dgym = dgy - dgtest

                update_status = update_trial_interval(stx, fxm, dgxm, sty, fym, &
                    dgym, step, fm, dgm, stmin, stmax, bracketed)
                fx = fxm + stx * dgtest
                fy = fym + sty * dgtest
                dgx = dgxm + dgtest
                dgy = dgym + dgtest
            else
                update_status = update_trial_interval(stx, fx, dgx, sty, fy, dgy, &
                    step, f, dg, stmin, stmax, bracketed)
            end if

            if (bracketed) then
                if (0.66_dp * previous_width <= abs(sty - stx)) then
                    step = stx + 0.5_dp * (sty - stx)
                end if
                previous_width = width
                width = abs(sty - stx)
            end if
        end do
    end function line_search_morethuente

    integer function update_trial_interval(x, fx, dx, y, fy, dy, t, ft, dt, &
            tmin, tmax, bracketed) result(status)
        real(dp), intent(inout) :: x, fx, dx, y, fy, dy, t, ft, dt
        real(dp), intent(in) :: tmin, tmax
        logical, intent(inout) :: bracketed

        logical :: derivative_sign_differs, bounded
        real(dp) :: cubic, quadratic, new_trial, limit

        status = lbfgs_success
        derivative_sign_differs = dt * dx < 0.0_dp

        if (bracketed) then
            if (t <= min(x, y) .or. t >= max(x, y)) then
                status = lbfgserr_outofinterval
                return
            else if (dx * (t - x) >= 0.0_dp) then
                status = lbfgserr_increasegradient
                return
            else if (tmax < tmin) then
                status = lbfgserr_incorrect_tminmax
                return
            end if
        end if

        if (fx < ft) then
            bracketed = .true.
            bounded = .true.
            cubic = cubic_minimizer(x, fx, dx, t, ft, dt)
            quadratic = quadratic_minimizer(x, fx, dx, t, ft)
            if (abs(cubic - x) < abs(quadratic - x)) then
                new_trial = cubic
            else
                new_trial = cubic + 0.5_dp * (quadratic - cubic)
            end if
        else if (derivative_sign_differs) then
            bracketed = .true.
            bounded = .false.
            cubic = cubic_minimizer(x, fx, dx, t, ft, dt)
            quadratic = secant_minimizer(x, dx, t, dt)
            if (abs(cubic - t) > abs(quadratic - t)) then
                new_trial = cubic
            else
                new_trial = quadratic
            end if
        else if (abs(dt) < abs(dx)) then
            bounded = .true.
            cubic = safeguarded_cubic_minimizer(x, fx, dx, t, ft, dt, tmin, tmax)
            quadratic = secant_minimizer(x, dx, t, dt)
            if (bracketed) then
                if (abs(t - cubic) < abs(t - quadratic)) then
                    new_trial = cubic
                else
                    new_trial = quadratic
                end if
            else
                if (abs(t - cubic) > abs(t - quadratic)) then
                    new_trial = cubic
                else
                    new_trial = quadratic
                end if
            end if
        else
            bounded = .false.
            if (bracketed) then
                new_trial = cubic_minimizer(t, ft, dt, y, fy, dy)
            else if (x < t) then
                new_trial = tmax
            else
                new_trial = tmin
            end if
        end if

        if (fx < ft) then
            y = t
            fy = ft
            dy = dt
        else
            if (derivative_sign_differs) then
                y = x
                fy = fx
                dy = dx
            end if
            x = t
            fx = ft
            dx = dt
        end if

        new_trial = max(tmin, min(tmax, new_trial))
        if (bracketed .and. bounded) then
            limit = x + 0.66_dp * (y - x)
            if (x < y) then
                new_trial = min(limit, new_trial)
            else
                new_trial = max(limit, new_trial)
            end if
        end if
        t = new_trial
    end function update_trial_interval

    pure real(dp) function cubic_minimizer(u, fu, du, v, fv, dv) result(value)
        real(dp), intent(in) :: u, fu, du, v, fv, dv
        real(dp) :: distance, theta, scale, gamma, ratio, p, q, radicand

        distance = v - u
        if (distance == 0.0_dp) then
            value = u
            return
        end if
        theta = 3.0_dp * (fu - fv) / distance + du + dv
        scale = max(abs(theta), abs(du), abs(dv))
        if (scale == 0.0_dp) then
            value = 0.5_dp * (u + v)
            return
        end if
        radicand = max(0.0_dp, (theta / scale)**2 - (du / scale) * (dv / scale))
        gamma = scale * sqrt(radicand)
        if (v < u) gamma = -gamma
        p = gamma - du + theta
        q = gamma - du + gamma + dv
        if (q == 0.0_dp) then
            value = 0.5_dp * (u + v)
        else
            ratio = p / q
            value = u + ratio * distance
        end if
    end function cubic_minimizer

    pure real(dp) function safeguarded_cubic_minimizer(u, fu, du, v, fv, dv, &
            xmin, xmax) result(value)
        real(dp), intent(in) :: u, fu, du, v, fv, dv, xmin, xmax
        real(dp) :: distance, theta, scale, gamma, ratio, p, q, a, radicand

        distance = v - u
        if (distance == 0.0_dp) then
            value = max(xmin, min(xmax, u))
            return
        end if
        theta = 3.0_dp * (fu - fv) / distance + du + dv
        scale = max(abs(theta), abs(du), abs(dv))
        if (scale == 0.0_dp) then
            value = max(xmin, min(xmax, 0.5_dp * (u + v)))
            return
        end if
        a = theta / scale
        radicand = max(0.0_dp, a * a - (du / scale) * (dv / scale))
        gamma = scale * sqrt(radicand)
        if (u < v) gamma = -gamma
        p = gamma - dv + theta
        q = gamma - dv + gamma + du
        if (q == 0.0_dp) then
            value = merge(xmax, xmin, a < 0.0_dp)
            return
        end if
        ratio = p / q
        if (ratio < 0.0_dp .and. gamma /= 0.0_dp) then
            value = v - ratio * distance
        else if (a < 0.0_dp) then
            value = xmax
        else
            value = xmin
        end if
    end function safeguarded_cubic_minimizer

    pure real(dp) function quadratic_minimizer(u, fu, du, v, fv) result(value)
        real(dp), intent(in) :: u, fu, du, v, fv
        real(dp) :: distance, denominator

        distance = v - u
        if (distance == 0.0_dp) then
            value = u
            return
        end if
        denominator = 2.0_dp * ((fu - fv) / distance + du)
        if (denominator == 0.0_dp) then
            value = 0.5_dp * (u + v)
        else
            value = u + du / denominator * distance
        end if
    end function quadratic_minimizer

    pure real(dp) function secant_minimizer(u, du, v, dv) result(value)
        real(dp), intent(in) :: u, du, v, dv
        real(dp) :: denominator

        denominator = dv - du
        if (denominator == 0.0_dp) then
            value = 0.5_dp * (u + v)
        else
            value = v + dv / denominator * (u - v)
        end if
    end function secant_minimizer

    pure real(dp) function owlqn_x1norm(x, start_index, end_index) result(norm)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: start_index, end_index

        norm = sum(abs(x(start_index:end_index)))
    end function owlqn_x1norm

    pure subroutine owlqn_pseudo_gradient(pg, x, g, coefficient, start_index, &
            end_index)
        real(dp), intent(out) :: pg(:)
        real(dp), intent(in) :: x(:), g(:), coefficient
        integer, intent(in) :: start_index, end_index
        integer :: i

        pg = g
        do i = start_index, end_index
            if (x(i) < 0.0_dp) then
                pg(i) = g(i) - coefficient
            else if (x(i) > 0.0_dp) then
                pg(i) = g(i) + coefficient
            else if (g(i) < -coefficient) then
                pg(i) = g(i) + coefficient
            else if (g(i) > coefficient) then
                pg(i) = g(i) - coefficient
            else
                pg(i) = 0.0_dp
            end if
        end do
    end subroutine owlqn_pseudo_gradient

    pure subroutine owlqn_project(x, sign_vector, start_index, end_index)
        real(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: sign_vector(:)
        integer, intent(in) :: start_index, end_index
        integer :: i

        do i = start_index, end_index
            if (x(i) * sign_vector(i) <= 0.0_dp) x(i) = 0.0_dp
        end do
    end subroutine owlqn_project

end module lbfgs_solver
