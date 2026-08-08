module test_max_problem
    use trust, only : dp
    implicit none
contains
    subroutine concave(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = 7.0_dp - (x(1)-2.0_dp)**2 - 3.0_dp*(x(2)+1.0_dp)**2
        gradient = [-2.0_dp*(x(1)-2.0_dp), -6.0_dp*(x(2)+1.0_dp)]
        hessian = 0.0_dp
        hessian(1,1) = -2.0_dp
        hessian(2,2) = -6.0_dp
        status = 0
    end subroutine concave
end module test_max_problem

program test_maximize
    use trust
    use test_max_problem, only : concave
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = [-4.0_dp, 5.0_dp]
    opt%rinit = 0.5_dp
    opt%rmax = 10.0_dp
    opt%minimize = .false.
    call trust_optimize(concave, x0, opt, res)
    if (.not. res%converged) error stop 'maximization did not converge'
    if (maxval(abs(res%argument - [2.0_dp, -1.0_dp])) > 1.0e-8_dp) error stop 'maximum wrong'
    if (abs(res%value - 7.0_dp) > 1.0e-10_dp) error stop 'maximum value wrong'
    print *, 'PASS test_maximize'
end program test_maximize
