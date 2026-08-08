module rosenbrock_example_problem
    use trust, only : dp
    implicit none
contains
    subroutine objective(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = 100.0_dp*(x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2
        gradient = [-400.0_dp*x(1)*(x(2)-x(1)**2)+2.0_dp*(x(1)-1.0_dp), &
            200.0_dp*(x(2)-x(1)**2)]
        hessian(1,1) = 1200.0_dp*x(1)**2 - 400.0_dp*x(2) + 2.0_dp
        hessian(1,2) = -400.0_dp*x(1)
        hessian(2,1) = hessian(1,2)
        hessian(2,2) = 200.0_dp
        status = 0
    end subroutine objective
end module rosenbrock_example_problem

program rosenbrock_example
    use trust
    use rosenbrock_example_problem, only : objective
    implicit none
    type(trust_options) :: options
    type(trust_result) :: result
    options%rinit = 1.0_dp
    options%rmax = 5.0_dp
    call trust_optimize(objective, [3.0_dp, 1.0_dp], options, result)
    print '(a,l1)', 'converged: ', result%converged
    print '(a,i0)', 'iterations: ', result%iterations
    print '(a,es16.8)', 'value: ', result%value
    print '(a,2f16.9)', 'argument: ', result%argument
end program rosenbrock_example
