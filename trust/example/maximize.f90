module maximize_example_problem
    use trust, only : dp
    implicit none
contains
    subroutine objective(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = 10.0_dp - (x(1)-2.0_dp)**2 - 4.0_dp*(x(2)+1.0_dp)**2
        gradient = [-2.0_dp*(x(1)-2.0_dp), -8.0_dp*(x(2)+1.0_dp)]
        hessian = 0.0_dp
        hessian(1,1) = -2.0_dp
        hessian(2,2) = -8.0_dp
        status = 0
    end subroutine objective
end module maximize_example_problem

program maximize_example
    use trust
    use maximize_example_problem, only : objective
    implicit none
    type(trust_options) :: options
    type(trust_result) :: result
    options%rinit = 0.5_dp
    options%rmax = 10.0_dp
    options%minimize = .false.
    call trust_optimize(objective, [-4.0_dp, 5.0_dp], options, result)
    print '(a,l1)', 'converged: ', result%converged
    print '(a,es16.8)', 'maximum: ', result%value
    print '(a,2f16.9)', 'argument: ', result%argument
end program maximize_example
