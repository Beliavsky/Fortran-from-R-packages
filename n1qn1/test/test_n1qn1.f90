module test_functions
    use n1qn1_module, only : dp
    implicit none
contains
    function banana_value(x, user_data) result(f)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: f
        integer :: i
        if (present(user_data)) continue
        f = 1.0_dp
        do i = 2, size(x)
            f = f + 100.0_dp * (x(i) - x(i - 1) ** 2) ** 2 + (1.0_dp - x(i)) ** 2
        end do
    end function banana_value

    subroutine banana_gradient(x, g, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(inout), optional :: user_data
        integer :: i, n
        if (present(user_data)) continue
        n = size(x)
        g = 0.0_dp
        g(1) = -400.0_dp * (x(2) - x(1) ** 2) * x(1)
        do i = 2, n - 1
            g(i) = 200.0_dp * (x(i) - x(i - 1) ** 2) &
                 - 400.0_dp * (x(i + 1) - x(i) ** 2) * x(i) &
                 - 2.0_dp * (1.0_dp - x(i))
        end do
        g(n) = 200.0_dp * (x(n) - x(n - 1) ** 2) - 2.0_dp * (1.0_dp - x(n))
    end subroutine banana_gradient

    function quadratic_value(x, user_data) result(f)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: f
        if (present(user_data)) continue
        f = 0.5_dp * (4.0_dp * (x(1) - 1.0_dp) ** 2 + &
            2.0_dp * (x(2) + 2.0_dp) ** 2 + 6.0_dp * (x(3) - 0.5_dp) ** 2)
    end function quadratic_value

    subroutine quadratic_gradient(x, g, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(inout), optional :: user_data
        if (present(user_data)) continue
        g = [4.0_dp * (x(1) - 1.0_dp), 2.0_dp * (x(2) + 2.0_dp), &
             6.0_dp * (x(3) - 0.5_dp)]
    end subroutine quadratic_gradient

    function square_value(x, user_data) result(f)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: f
        if (present(user_data)) continue
        f = x(1) ** 2
    end function square_value

    subroutine square_gradient(x, g, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(inout), optional :: user_data
        if (present(user_data)) continue
        g(1) = 2.0_dp * x(1)
    end subroutine square_gradient

    subroutine stop_immediately(iteration, evaluations, x, value, gradient, stop, user_data)
        integer, intent(in) :: iteration, evaluations
        real(dp), intent(in) :: x(:), value, gradient(:)
        logical, intent(out) :: stop
        class(*), intent(inout), optional :: user_data
        stop = iteration >= 1 .and. evaluations >= 1 .and. size(x) == size(gradient) &
            .and. value < huge(1.0_dp)
        if (present(user_data)) then
            select type (user_data)
            type is (integer)
                user_data = 1
            end select
        end if
    end subroutine stop_immediately
end module test_functions

program test_n1qn1
    use n1qn1_module, only : dp, n1qn1_control_t, n1qn1_result_t, n1qn1_minimize, &
                             n1qn1_success, n1qn1_max_evaluations
    use test_functions
    implicit none

    integer :: failures
    failures = 0
    call test_banana(failures)
    call test_initial_hessian(failures)
    call test_reuse_curvature(failures)
    call test_quadratic(failures)
    call test_single_variable(failures)
    call test_limits(failures)
    call test_progress_stop(failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'FAILED tests: ', failures
        error stop 1
    end if
    write(*, '(a)') 'All n1qn1 tests passed.'

contains
    subroutine test_banana(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: result
        type(n1qn1_control_t) :: control
        real(dp) :: x0(3)

        x0 = [1.02_dp, 1.02_dp, 1.02_dp]
        control%epsilon = epsilon(1.0_dp)
        control%max_iterations = 100
        control%max_evaluations = 100
        call n1qn1_minimize(banana_value, banana_gradient, x0, result, control)
        call check(maxval(abs(result%x - 1.0_dp)) < 1.0e-10_dp, 'banana x', failures)
        call check(abs(result%value - 1.0_dp) < 1.0e-12_dp, 'banana value', failures)
        call check(result%function_evaluations == 40, 'banana evaluation count', failures)
        call check(maxval(abs(result%hessian - reshape([ &
            799.995909385953_dp, -399.614075225545_dp, -0.196213499712260_dp, &
            -399.614075225545_dp, 1002.594973913260_dp, -400.319316558786_dp, &
            -0.196213499712260_dp, -400.319316558786_dp, 202.170906236328_dp], [3, 3]))) &
            < 2.0e-9_dp, 'banana Hessian', failures)
    end subroutine test_banana

    subroutine test_initial_hessian(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: result
        real(dp) :: x0(3), h0(3, 3)

        x0 = [1.02_dp, 1.02_dp, 1.02_dp]
        h0 = reshape([ &
            797.861115_dp, -393.801473_dp, -2.795134_dp, &
            -393.801473_dp, 991.271179_dp, -395.382900_dp, &
            -2.795134_dp, -395.382900_dp, 200.024349_dp], [3, 3])
        call n1qn1_minimize(banana_value, banana_gradient, x0, result, initial_hessian=h0)
        call check(maxval(abs(result%x - 1.0_dp)) < 1.0e-9_dp, 'initial Hessian x', failures)
        call check(result%function_evaluations == 29, 'initial Hessian evaluations', failures)
        call check(maxval(abs(result%hessian - reshape([ &
            800.030807707827_dp, -399.878160447993_dp, -0.0526692400192971_dp, &
            -399.878160447993_dp, 1001.84045503084_dp, -399.890537542132_dp, &
            -0.0526692400192971_dp, -399.890537542132_dp, 201.932617669621_dp], [3, 3]))) &
            < 2.0e-9_dp, 'initial Hessian output', failures)
    end subroutine test_initial_hessian

    subroutine test_reuse_curvature(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: first, second, restarted
        real(dp) :: x0(3)

        x0 = [1.02_dp, 1.02_dp, 1.02_dp]
        call n1qn1_minimize(banana_value, banana_gradient, x0, first)
        call n1qn1_minimize(banana_value, banana_gradient, x0, second, &
                            initial_hessian=first%hessian)
        call n1qn1_minimize(banana_value, banana_gradient, x0, restarted, &
                            initial_factor=first%factor)
        call check(second%function_evaluations == 33, 'reused full Hessian evaluations', failures)
        call check(maxval(abs(second%x - 1.0_dp)) < 1.0e-10_dp, 'reused full Hessian x', failures)
        call check(maxval(abs(restarted%x - 1.0_dp)) < 1.0e-10_dp, 'factor restart x', failures)
        call check(restarted%function_evaluations < first%function_evaluations, &
                   'factor restart reduces evaluations', failures)
    end subroutine test_reuse_curvature

    subroutine test_quadratic(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: result
        real(dp) :: x0(3)
        x0 = [4.0_dp, -5.0_dp, 2.0_dp]
        call n1qn1_minimize(quadratic_value, quadratic_gradient, x0, result)
        call check(maxval(abs(result%x - [1.0_dp, -2.0_dp, 0.5_dp])) < 1.0e-10_dp, &
                   'quadratic solution', failures)
        call check(result%status == n1qn1_success, 'quadratic status', failures)
    end subroutine test_quadratic

    subroutine test_single_variable(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: result
        call n1qn1_minimize(square_value, square_gradient, [3.0_dp], result)
        call check(abs(result%x(1)) < 1.0e-12_dp, 'single variable x', failures)
        call check(abs(result%value) < 1.0e-20_dp, 'single variable value', failures)
    end subroutine test_single_variable

    subroutine test_limits(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: result
        type(n1qn1_control_t) :: control
        control%max_evaluations = 2
        call n1qn1_minimize(banana_value, banana_gradient, [1.02_dp, 1.02_dp, 1.02_dp], &
                            result, control)
        call check(result%status == n1qn1_max_evaluations, 'evaluation limit', failures)
    end subroutine test_limits

    subroutine test_progress_stop(failures)
        integer, intent(inout) :: failures
        type(n1qn1_result_t) :: result
        integer :: marker
        marker = 0
        call n1qn1_minimize(quadratic_value, quadratic_gradient, [4.0_dp, -5.0_dp, 2.0_dp], &
                            result, user_data=marker, progress=stop_immediately)
        call check(marker == 1, 'progress callback data', failures)
        call check(result%function_evaluations == 1, 'progress stop evaluations', failures)
    end subroutine test_progress_stop

    subroutine check(condition, label, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        if (.not. condition) then
            failures = failures + 1
            write(*, '(a)') 'FAIL: ' // trim(label)
        end if
    end subroutine check
end program test_n1qn1
