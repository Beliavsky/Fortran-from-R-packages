program test_nonneg_cg
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use nonneg_cg
    implicit none

    type :: shift_data_t
        real(dp) :: target
    end type shift_data_t

    integer :: failures

    failures = 0
    call test_rosenbrock_reference(failures)
    call test_rosenbrock_convergence(failures)
    call test_active_boundary(failures)
    call test_boundary_stationary(failures)
    call test_user_data(failures)
    call test_max_iteration(failures)
    call test_max_evaluation(failures)
    call test_cancellation(failures)
    call test_invalid_input(failures)
    call test_nonfinite_gradient(failures)

    if (failures /= 0) then
        write (*, '(a,i0)') 'FAILED tests: ', failures
        error stop 1
    end if
    write (*, '(a)') 'All nonneg-cg tests passed.'

contains

    subroutine check(condition, name, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name
        integer, intent(inout) :: failures

        if (.not. condition) then
            failures = failures + 1
            write (*, '(a,a)') 'FAIL: ', name
        end if
    end subroutine check

    subroutine test_rosenbrock_reference(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_control_t) :: control
        type(nonneg_cg_result_t) :: result

        x = [0.0_dp, 2.0_dp]
        control%tol = 1.0e-8_dp
        control%maxiter = 5
        call minimize_nonneg_cg(x, rosenbrock, rosenbrock_gradient, result, control)

        call check(result%status == nonneg_cg_stop_maxiter, 'C reference status', failures)
        call check(result%niter == 5, 'C reference iteration count', failures)
        call check(result%nfeval == 24, 'C reference legacy evaluation count', failures)
        call check(result%objective_calls == 29, 'C reference true objective calls', failures)
        call check(maxval(abs(result%x - [1.1822794465797628_dp, &
            1.4201196010889658_dp])) < 5.0e-13_dp, 'C reference parameters', failures)
        call check(abs(result%fun - 0.083110622851909594_dp) < 5.0e-13_dp, &
            'C reference objective', failures)
    end subroutine test_rosenbrock_reference

    subroutine test_rosenbrock_convergence(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_control_t) :: control
        type(nonneg_cg_result_t) :: result

        x = [0.0_dp, 2.0_dp]
        control%tol = 1.0e-8_dp
        call minimize_nonneg_cg(x, rosenbrock, rosenbrock_gradient, result, control)

        call check(result%status == nonneg_cg_tol_achieved, 'Rosenbrock convergence status', failures)
        call check(maxval(abs(result%x - 1.0_dp)) < 1.0e-3_dp, &
            'Rosenbrock converged parameters', failures)
        call check(result%fun < 1.0e-8_dp, 'Rosenbrock converged objective', failures)
        call check(result%objective_calls == result%nfeval + result%niter, &
            'objective-call accounting', failures)
    end subroutine test_rosenbrock_convergence

    subroutine test_active_boundary(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_control_t) :: control
        type(nonneg_cg_result_t) :: result

        x = [2.0_dp, 0.0_dp]
        control%tol = 1.0e-12_dp
        control%maxiter = 500
        call minimize_nonneg_cg(x, boundary_objective, boundary_gradient, result, control)

        call check(result%status == nonneg_cg_tol_achieved, 'boundary status', failures)
        call check(abs(result%x(1)) < 1.0e-12_dp, 'active lower bound', failures)
        call check(abs(result%x(2) - 3.0_dp) < 1.0e-8_dp, 'free variable optimum', failures)
        call check(all(result%x >= 0.0_dp), 'boundary feasibility', failures)
    end subroutine test_active_boundary

    subroutine test_boundary_stationary(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_result_t) :: result

        x = 0.0_dp
        call minimize_nonneg_cg(x, positive_linear, positive_linear_gradient, result)
        call check(result%status == nonneg_cg_tol_achieved, 'stationary boundary status', failures)
        call check(result%niter == 0, 'stationary boundary iterations', failures)
        call check(maxval(abs(result%x)) <= tiny(1.0_dp), 'stationary boundary point', failures)
    end subroutine test_boundary_stationary

    subroutine test_user_data(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(1)
        type(shift_data_t) :: data
        type(nonneg_cg_control_t) :: control
        type(nonneg_cg_result_t) :: result

        data%target = 4.0_dp
        x = [1.0_dp]
        control%tol = 1.0e-12_dp
        call minimize_nonneg_cg(x, shifted_objective, shifted_gradient, result, control, data)
        call check(abs(result%x(1) - 4.0_dp) < 1.0e-10_dp, 'user data target', failures)
    end subroutine test_user_data

    subroutine test_max_iteration(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_control_t) :: control
        type(nonneg_cg_result_t) :: result

        x = [0.0_dp, 2.0_dp]
        control%maxiter = 1
        control%tol = 0.0_dp
        call minimize_nonneg_cg(x, rosenbrock, rosenbrock_gradient, result, control)
        call check(result%status == nonneg_cg_stop_maxiter, 'maximum iteration status', failures)
        call check(result%niter == 1, 'maximum iteration count', failures)
    end subroutine test_max_iteration

    subroutine test_max_evaluation(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_control_t) :: control
        type(nonneg_cg_result_t) :: result

        x = [0.0_dp, 2.0_dp]
        control%maxnfeval = 2
        control%tol = 0.0_dp
        call minimize_nonneg_cg(x, rosenbrock, rosenbrock_gradient, result, control)
        call check(result%status == nonneg_cg_stop_maxnfeval, 'maximum evaluation status', failures)
        call check(maxval(abs(result%x - [0.01_dp, 0.0_dp])) < 2.0e-15_dp, &
            'maximum evaluation reversion', failures)
    end subroutine test_max_evaluation

    subroutine test_cancellation(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(2)
        type(nonneg_cg_result_t) :: result

        x = [0.0_dp, 2.0_dp]
        call minimize_nonneg_cg(x, rosenbrock, rosenbrock_gradient, result, monitor=cancel_after_one)
        call check(result%status == nonneg_cg_cancelled, 'monitor cancellation status', failures)
        call check(result%niter == 1, 'monitor cancellation iteration', failures)
    end subroutine test_cancellation

    subroutine test_invalid_input(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(1)
        type(nonneg_cg_result_t) :: result

        x = [-1.0_dp]
        call minimize_nonneg_cg(x, positive_linear, positive_linear_gradient, result)
        call check(result%status == nonneg_cg_invalid_input, 'invalid starting point', failures)
    end subroutine test_invalid_input

    subroutine test_nonfinite_gradient(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(1)
        type(nonneg_cg_result_t) :: result

        x = [1.0_dp]
        call minimize_nonneg_cg(x, positive_linear, nan_gradient, result)
        call check(result%status == nonneg_cg_nonfinite_value, 'nonfinite gradient status', failures)
    end subroutine test_nonfinite_gradient

    function rosenbrock(x, user_data) result(value)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: value

        value = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
    end function rosenbrock

    subroutine rosenbrock_gradient(x, gradient, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: gradient(:)
        class(*), intent(inout), optional :: user_data

        gradient(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - 2.0_dp * (1.0_dp - x(1))
        gradient(2) = 200.0_dp * (x(2) - x(1)**2)
    end subroutine rosenbrock_gradient

    function boundary_objective(x, user_data) result(value)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: value

        value = 0.5_dp * ((x(1) + 1.0_dp)**2 + (x(2) - 3.0_dp)**2)
    end function boundary_objective

    subroutine boundary_gradient(x, gradient, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: gradient(:)
        class(*), intent(inout), optional :: user_data

        gradient = [x(1) + 1.0_dp, x(2) - 3.0_dp]
    end subroutine boundary_gradient

    function positive_linear(x, user_data) result(value)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: value

        value = sum(x)
    end function positive_linear

    subroutine positive_linear_gradient(x, gradient, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: gradient(:)
        class(*), intent(inout), optional :: user_data

        gradient = 1.0_dp
    end subroutine positive_linear_gradient

    function shifted_objective(x, user_data) result(value)
        real(dp), intent(in) :: x(:)
        class(*), intent(inout), optional :: user_data
        real(dp) :: value

        select type (data => user_data)
        type is (shift_data_t)
            value = 0.5_dp * (x(1) - data%target)**2
        class default
            value = huge(1.0_dp)
        end select
    end function shifted_objective

    subroutine shifted_gradient(x, gradient, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: gradient(:)
        class(*), intent(inout), optional :: user_data

        select type (data => user_data)
        type is (shift_data_t)
            gradient(1) = x(1) - data%target
        class default
            gradient = huge(1.0_dp)
        end select
    end subroutine shifted_gradient

    function cancel_after_one(x, value, iteration, user_data) result(cancel)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: value
        integer, intent(in) :: iteration
        class(*), intent(inout), optional :: user_data
        logical :: cancel

        cancel = iteration >= 1
    end function cancel_after_one

    subroutine nan_gradient(x, gradient, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: gradient(:)
        class(*), intent(inout), optional :: user_data

        gradient = ieee_value(0.0_dp, ieee_quiet_nan)
    end subroutine nan_gradient

end program test_nonneg_cg
