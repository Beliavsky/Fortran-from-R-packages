module test_scaling_problem
    use trust, only : dp
    implicit none
contains
    subroutine quadratic(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = (x(1)-2.0_dp)**2 + 25.0_dp*(x(2)+3.0_dp)**2
        gradient = [2.0_dp*(x(1)-2.0_dp), 50.0_dp*(x(2)+3.0_dp)]
        hessian = 0.0_dp
        hessian(1,1) = 2.0_dp
        hessian(2,2) = 50.0_dp
        status = 0
    end subroutine quadratic
end module test_scaling_problem

program test_scaling
    use trust
    use test_scaling_problem, only : quadratic
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = [10.0_dp, 10.0_dp]
    opt%rinit = 1.0_dp
    opt%rmax = 100.0_dp
    allocate(opt%parscale(2))
    opt%parscale = [5.0_dp, 0.5_dp]
    call trust_optimize(quadratic, x0, opt, res)
    if (.not. res%converged) error stop 'scaled problem did not converge'
    if (maxval(abs(res%argument - [2.0_dp, -3.0_dp])) > 1.0e-8_dp) error stop 'scaled solution wrong'
    print *, 'PASS test_scaling'
end program test_scaling
