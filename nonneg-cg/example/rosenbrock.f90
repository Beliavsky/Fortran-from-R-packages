program rosenbrock_example
    use nonneg_cg, only : dp, minimize_nonneg_cg, nonneg_cg_control_t, nonneg_cg_result_t
    implicit none

    real(dp) :: x(2)
    type(nonneg_cg_control_t) :: control
    type(nonneg_cg_result_t) :: result

    x = [0.0_dp, 2.0_dp]
    control%tol = 1.0e-8_dp
    control%verbose = .true.

    call minimize_nonneg_cg(x, rosenbrock, rosenbrock_gradient, result, control)

    write (*, '(a,2f16.10)') 'x = ', result%x
    write (*, '(a,es16.8)') 'f = ', result%fun
    write (*, '(a,i0)') 'iterations = ', result%niter
    write (*, '(a,i0)') 'legacy nfeval = ', result%nfeval
    write (*, '(a,i0)') 'actual objective calls = ', result%objective_calls

contains

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

end program rosenbrock_example
