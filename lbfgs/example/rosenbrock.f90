program rosenbrock_example
    use lbfgs
    implicit none

    real(dp) :: x(2)
    type(lbfgs_parameter_t) :: parameters
    type(lbfgs_result_t) :: result

    x = [-1.2_dp, 1.0_dp]
    parameters = lbfgs_parameter_t()
    parameters%epsilon = 1.0e-9_dp

    call lbfgs_minimize(rosenbrock, x, result, parameters)

    print '(a,2(1x,f14.9))', 'solution:', x
    print '(a,1x,es16.8)', 'objective:', result%value
    print '(a,1x,i0)', 'iterations:', result%iterations
    print '(a,1x,i0)', 'evaluations:', result%evaluations
    print '(a,1x,a)', 'status:', result%message

contains

    subroutine rosenbrock(x, f, g, step, user_data)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(:)
        real(dp), intent(in) :: step
        class(*), intent(inout), optional :: user_data

        f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
        g(1) = -400.0_dp * x(1) * (x(2) - x(1)**2) - &
            2.0_dp * (1.0_dp - x(1))
        g(2) = 200.0_dp * (x(2) - x(1)**2)
    end subroutine rosenbrock

end program rosenbrock_example
