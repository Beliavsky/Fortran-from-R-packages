program test_lbfgs
    use lbfgs
    implicit none

    type :: quadratic_data_t
        real(dp) :: target(3)
    end type quadratic_data_t

    integer :: failures

    failures = 0
    call test_quadratic_linesearches(failures)
    call test_rosenbrock(failures)
    call test_owlqn(failures)
    call test_partial_owlqn(failures)
    call test_separate_callbacks(failures)
    call test_user_data(failures)
    call test_progress_cancel(failures)
    call test_validation(failures)
    call test_max_iterations(failures)
    call test_already_minimized(failures)

    if (failures /= 0) then
        print '(a,i0)', 'FAILED tests: ', failures
        error stop 1
    end if
    print '(a)', 'All lbfgs tests passed.'

contains

    subroutine test_quadratic_linesearches(failures)
        integer, intent(inout) :: failures
        integer, parameter :: algorithms(4) = [ &
            lbfgs_linesearch_morethuente, &
            lbfgs_linesearch_backtracking_armijo, &
            lbfgs_linesearch_backtracking_wolfe, &
            lbfgs_linesearch_backtracking_strong_wolfe]
        real(dp) :: x(20), target(20)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result
        integer :: i

        target = [(0.1_dp * real(i, dp), i = 1, 20)]
        do i = 1, size(algorithms)
            x = -2.0_dp
            param = lbfgs_parameter_t()
            param%linesearch = algorithms(i)
            param%epsilon = 1.0e-10_dp
            param%max_iterations = 300
            call lbfgs_minimize(quadratic_evaluate, x, result, param)
            call assert_true(result%status == lbfgs_success, &
                'quadratic line-search status', failures)
            call assert_close(maxval(abs(x - target)), 0.0_dp, 1.0e-7_dp, &
                'quadratic solution', failures)
            call assert_true(result%value < 1.0e-12_dp, &
                'quadratic objective', failures)
        end do
    end subroutine test_quadratic_linesearches

    subroutine test_rosenbrock(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        x = [-1.2_dp, 1.0_dp]
        param = lbfgs_parameter_t()
        param%epsilon = 1.0e-9_dp
        param%max_iterations = 1000
        call lbfgs_minimize(rosenbrock_evaluate, x, result, param)
        call assert_true(result%status == lbfgs_success, &
            'Rosenbrock status', failures)
        call assert_close(maxval(abs(x - 1.0_dp)), 0.0_dp, 2.0e-5_dp, &
            'Rosenbrock solution', failures)
        call assert_true(result%value < 1.0e-12_dp, &
            'Rosenbrock objective', failures)
    end subroutine test_rosenbrock

    subroutine test_owlqn(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(4), expected(4)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        x = 0.0_dp
        expected = [2.5_dp, -0.5_dp, 0.0_dp, -3.5_dp]
        param = lbfgs_parameter_t()
        param%orthantwise_c = 0.5_dp
        param%linesearch = lbfgs_linesearch_backtracking
        param%epsilon = 1.0e-10_dp
        param%max_iterations = 200
        call lbfgs_minimize(lasso_evaluate, x, result, param)
        call assert_true(result%status == lbfgs_success, &
            'OWL-QN status', failures)
        call assert_close(maxval(abs(x - expected)), 0.0_dp, 1.0e-7_dp, &
            'OWL-QN soft threshold', failures)
    end subroutine test_owlqn

    subroutine test_partial_owlqn(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(4), expected(4)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        x = 0.0_dp
        expected = [3.0_dp, -0.5_dp, 0.0_dp, -4.0_dp]
        param = lbfgs_parameter_t()
        param%orthantwise_c = 0.5_dp
        param%orthantwise_start = 2
        param%orthantwise_end = 3
        param%linesearch = lbfgs_linesearch_backtracking
        param%epsilon = 1.0e-10_dp
        param%max_iterations = 200
        call lbfgs_minimize(lasso_evaluate, x, result, param)
        call assert_true(result%status == lbfgs_success, &
            'partial OWL-QN status', failures)
        call assert_close(maxval(abs(x - expected)), 0.0_dp, 1.0e-7_dp, &
            'partial OWL-QN solution', failures)
    end subroutine test_partial_owlqn

    subroutine test_separate_callbacks(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(20), target(20)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result
        integer :: i

        target = [(0.1_dp * real(i, dp), i = 1, 20)]
        x = 3.0_dp
        param = lbfgs_parameter_t()
        param%epsilon = 1.0e-10_dp
        call lbfgs_minimize(quadratic_objective, quadratic_gradient, x, result, param)
        call assert_true(result%status == lbfgs_success, &
            'separate callbacks status', failures)
        call assert_close(maxval(abs(x - target)), 0.0_dp, 1.0e-7_dp, &
            'separate callbacks solution', failures)
    end subroutine test_separate_callbacks


    subroutine test_user_data(failures)
        integer, intent(inout) :: failures
        type(quadratic_data_t) :: data
        real(dp) :: x(3)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        data%target = [1.5_dp, -2.0_dp, 0.25_dp]
        x = 0.0_dp
        param = lbfgs_parameter_t()
        param%epsilon = 1.0e-11_dp
        call lbfgs_minimize(data_evaluate, x, result, param, user_data=data)
        call assert_true(result%status == lbfgs_success, &
            'user-data callback status', failures)
        call assert_close(maxval(abs(x - data%target)), 0.0_dp, 1.0e-8_dp, &
            'user-data callback solution', failures)
    end subroutine test_user_data

    subroutine test_progress_cancel(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        x = [-1.2_dp, 1.0_dp]
        param = lbfgs_parameter_t()
        param%max_iterations = 1000
        call lbfgs_minimize(rosenbrock_evaluate, x, result, param, cancel_progress)
        call assert_true(result%status == lbfgserr_canceled, &
            'progress cancellation status', failures)
        call assert_true(result%iterations <= 4, &
            'progress cancellation iteration', failures)
    end subroutine test_progress_cancel

    subroutine test_validation(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        x = 0.0_dp
        param = lbfgs_parameter_t()
        param%orthantwise_c = 1.0_dp
        param%linesearch = lbfgs_linesearch_morethuente
        call lbfgs_minimize(rosenbrock_evaluate, x, result, param)
        call assert_true(result%status == lbfgserr_invalid_linesearch, &
            'invalid OWL-QN line search', failures)

        param = lbfgs_parameter_t()
        param%orthantwise_start = 3
        call lbfgs_minimize(rosenbrock_evaluate, x, result, param)
        call assert_true(result%status == lbfgserr_invalid_orthantwise_start, &
            'invalid OWL-QN start', failures)
    end subroutine test_validation

    subroutine test_max_iterations(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result

        x = [-1.2_dp, 1.0_dp]
        param = lbfgs_parameter_t()
        param%max_iterations = 1
        call lbfgs_minimize(rosenbrock_evaluate, x, result, param)
        call assert_true(result%status == lbfgserr_maximumiteration, &
            'maximum iteration status', failures)
    end subroutine test_max_iterations

    subroutine test_already_minimized(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(20)
        type(lbfgs_parameter_t) :: param
        type(lbfgs_result_t) :: result
        integer :: i

        x = [(0.1_dp * real(i, dp), i = 1, 20)]
        param = lbfgs_parameter_t()
        call lbfgs_minimize(quadratic_evaluate, x, result, param)
        call assert_true(result%status == lbfgs_already_minimized, &
            'already minimized status', failures)
    end subroutine test_already_minimized

    subroutine quadratic_evaluate(x, f, g, step, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(:)
        real(dp), intent(in) :: step
        class(*), intent(inout), optional :: user_data
        real(dp) :: target(size(x))
        integer :: i

        target = [(0.1_dp * real(i, dp), i = 1, size(x))]
        g = x - target
        f = 0.5_dp * dot_product(g, g)
    end subroutine quadratic_evaluate

    subroutine quadratic_objective(x, f, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        class(*), intent(inout), optional :: user_data
        real(dp) :: target(size(x))
        integer :: i

        target = [(0.1_dp * real(i, dp), i = 1, size(x))]
        f = 0.5_dp * sum((x - target)**2)
    end subroutine quadratic_objective

    subroutine quadratic_gradient(x, g, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: target(size(x))
        integer :: i

        target = [(0.1_dp * real(i, dp), i = 1, size(x))]
        g = x - target
    end subroutine quadratic_gradient

    subroutine rosenbrock_evaluate(x, f, g, step, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(:)
        real(dp), intent(in) :: step
        class(*), intent(inout), optional :: user_data

        f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
        g(1) = -2.0_dp * (1.0_dp - x(1)) - &
            400.0_dp * x(1) * (x(2) - x(1)**2)
        g(2) = 200.0_dp * (x(2) - x(1)**2)
    end subroutine rosenbrock_evaluate

    subroutine lasso_evaluate(x, f, g, step, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(:)
        real(dp), intent(in) :: step
        class(*), intent(inout), optional :: user_data
        real(dp), parameter :: target(4) = [3.0_dp, -1.0_dp, 0.2_dp, -4.0_dp]

        g = x - target
        f = 0.5_dp * dot_product(g, g)
    end subroutine lasso_evaluate

    integer function cancel_progress(x, g, f, xnorm, gnorm, step, iteration, &
            line_evaluations, user_data) result(cancel)
        real(dp), intent(in) :: x(:), g(:), f, xnorm, gnorm, step
        integer, intent(in) :: iteration, line_evaluations
        class(*), intent(inout), optional :: user_data

        cancel = merge(1, 0, iteration >= 3)
    end function cancel_progress

    subroutine data_evaluate(x_current, f, g, step, user_data)
        real(dp), intent(in) :: x_current(:), step
        real(dp), intent(out) :: f, g(:)
        class(*), intent(inout), optional :: user_data

        select type (user_data)
        type is (quadratic_data_t)
            g = x_current - user_data%target
            f = 0.5_dp * dot_product(g, g)
        class default
            f = huge(1.0_dp)
            g = 0.0_dp
        end select
    end subroutine data_evaluate

    subroutine assert_true(condition, label, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        if (.not. condition) then
            print '(a)', 'FAIL: ' // trim(label)
            failures = failures + 1
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, label, failures)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures

        if (abs(actual - expected) > tolerance) then
            print '(a,2(1x,es16.8))', 'FAIL: ' // trim(label), actual, expected
            failures = failures + 1
        end if
    end subroutine assert_close

end program test_lbfgs
