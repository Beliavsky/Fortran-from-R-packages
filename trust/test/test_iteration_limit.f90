module test_unbounded_problem
    use trust, only : dp
    implicit none
contains
    subroutine saddle(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = x(1)**2 - x(2)**2
        gradient = [2.0_dp*x(1), -2.0_dp*x(2)]
        hessian = 0.0_dp
        hessian(1,1) = 2.0_dp
        hessian(2,2) = -2.0_dp
        status = 0
    end subroutine saddle
end module test_unbounded_problem

program test_iteration_limit
    use trust
    use test_unbounded_problem, only : saddle
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = [3.0_dp, 1.0_dp]
    opt%rinit = 1.0_dp
    opt%rmax = 5.0_dp
    opt%iterlim = 20
    call trust_optimize(saddle, x0, opt, res)
    if (res%converged) error stop 'unbounded problem should not converge'
    if (res%iterations /= 20) error stop 'iteration limit not respected'
    print *, 'PASS test_iteration_limit'
end program test_iteration_limit
