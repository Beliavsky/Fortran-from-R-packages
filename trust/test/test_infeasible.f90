module test_infeasible_problem
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use trust, only : dp
    implicit none
contains
    subroutine barrier(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: q, omq
        real(dp), parameter :: mu(2) = [1.0_dp, 2.0_dp]
        q = sum(x*x)
        status = 0
        if (q >= 1.0_dp) then
            value = ieee_value(0.0_dp, ieee_positive_inf)
            gradient = 0.0_dp
            hessian = 0.0_dp
            return
        end if
        omq = 1.0_dp - q
        value = dot_product(x, mu) - log(omq)
        gradient = mu + 2.0_dp*x/omq
        hessian = 4.0_dp*spread(x,2,2)*spread(x,1,2)/(omq*omq)
        hessian(1,1) = hessian(1,1) + 2.0_dp/omq
        hessian(2,2) = hessian(2,2) + 2.0_dp/omq
    end subroutine barrier
end module test_infeasible_problem

program test_infeasible
    use trust
    use test_infeasible_problem, only : barrier
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = 0.0_dp
    opt%rinit = 1.0_dp
    opt%rmax = 100.0_dp
    call trust_optimize(barrier, x0, opt, res)
    if (.not. res%converged) error stop 'barrier problem did not converge'
    if (sum(res%argument**2) >= 1.0_dp) error stop 'barrier solution infeasible'
    if (maxval(abs(res%gradient)) > 2.0e-6_dp) error stop 'barrier gradient too large'
    print *, 'PASS test_infeasible'
end program test_infeasible
