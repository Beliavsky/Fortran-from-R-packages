module test_hard_problem
    use trust, only : dp
    implicit none
contains
    subroutine quartic_saddle(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = x(1)**2 - x(2)**2 + 0.01_dp * (x(1)**4 + x(2)**4)
        gradient = [2.0_dp*x(1) + 0.04_dp*x(1)**3, -2.0_dp*x(2) + 0.04_dp*x(2)**3]
        hessian = 0.0_dp
        hessian(1,1) = 2.0_dp + 0.12_dp*x(1)**2
        hessian(2,2) = -2.0_dp + 0.12_dp*x(2)**2
        status = 0
    end subroutine quartic_saddle
end module test_hard_problem

program test_hard_case
    use trust
    use test_hard_problem, only : quartic_saddle
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = 0.0_dp
    opt%rinit = 1.0_dp
    opt%rmax = 5.0_dp
    opt%save_history = .true.
    call trust_optimize(quartic_saddle, x0, opt, res)
    if (.not. res%converged) error stop 'hard case did not converge'
    if (abs(res%argument(1)) > 1.0e-6_dp) error stop 'hard-case x1 wrong'
    if (abs(abs(res%argument(2)) - sqrt(50.0_dp)) > 2.0e-5_dp) error stop 'hard-case x2 wrong'
    if (abs(res%value + 25.0_dp) > 2.0e-8_dp) error stop 'hard-case value wrong'
    if (.not. any(res%history%step_type(:res%history%n) == trust_step_hard_hard)) &
        error stop 'hard-hard step not exercised'
    print *, 'PASS test_hard_case'
end program test_hard_case
