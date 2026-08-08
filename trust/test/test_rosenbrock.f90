module test_rosenbrock_problem
    use trust, only : dp
    implicit none
contains
    subroutine rosenbrock(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = 100.0_dp * (x(2) - x(1) * x(1)) ** 2 + (1.0_dp - x(1)) ** 2
        gradient(1) = -400.0_dp * x(1) * (x(2) - x(1) * x(1)) + 2.0_dp * (x(1) - 1.0_dp)
        gradient(2) = 200.0_dp * (x(2) - x(1) * x(1))
        hessian(1, 1) = 1200.0_dp * x(1) * x(1) - 400.0_dp * x(2) + 2.0_dp
        hessian(1, 2) = -400.0_dp * x(1)
        hessian(2, 1) = hessian(1, 2)
        hessian(2, 2) = 200.0_dp
        status = 0
    end subroutine rosenbrock
end module test_rosenbrock_problem

program test_rosenbrock
    use trust
    use test_rosenbrock_problem, only : rosenbrock
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = [3.0_dp, 1.0_dp]
    opt%rinit = 1.0_dp
    opt%rmax = 5.0_dp
    opt%save_history = .true.
    call trust_optimize(rosenbrock, x0, opt, res)
    if (.not. res%converged) error stop 'Rosenbrock did not converge'
    if (maxval(abs(res%argument - 1.0_dp)) > 2.0e-6_dp) error stop 'wrong Rosenbrock solution'
    if (res%value > 1.0e-10_dp) error stop 'wrong Rosenbrock value'
    if (res%history%n /= res%iterations) error stop 'history length mismatch'
    print *, 'PASS test_rosenbrock'
end program test_rosenbrock
