module test_step_cases_problem
    use trust, only : dp
    implicit none
contains
    subroutine easy_easy_obj(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = -x(1)**2 + x(2)**2 + x(1)
        gradient = [-2.0_dp*x(1) + 1.0_dp, 2.0_dp*x(2)]
        hessian = 0.0_dp
        hessian(1,1) = -2.0_dp
        hessian(2,2) = 2.0_dp
        status = 0
    end subroutine easy_easy_obj

    subroutine hard_easy_obj(x, value, gradient, hessian, status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: value, gradient(:), hessian(:, :)
        integer, intent(out) :: status
        value = -x(1)**2 + x(2)**2 + 10.0_dp*x(2)
        gradient = [-2.0_dp*x(1), 2.0_dp*x(2) + 10.0_dp]
        hessian = 0.0_dp
        hessian(1,1) = -2.0_dp
        hessian(2,2) = 2.0_dp
        status = 0
    end subroutine hard_easy_obj
end module test_step_cases_problem

program test_step_cases
    use trust
    use test_step_cases_problem, only : easy_easy_obj, hard_easy_obj
    implicit none
    type(trust_options) :: opt
    type(trust_result) :: res
    real(dp) :: x0(2)
    x0 = 0.0_dp
    opt%rinit = 1.0_dp
    opt%rmax = 5.0_dp
    opt%iterlim = 1
    opt%save_history = .true.

    call trust_optimize(easy_easy_obj, x0, opt, res)
    if (res%history%step_type(1) /= trust_step_easy_easy) error stop 'easy-easy not selected'

    call trust_optimize(hard_easy_obj, x0, opt, res)
    if (res%history%step_type(1) /= trust_step_hard_easy) error stop 'hard-easy not selected'

    print *, 'PASS test_step_cases'
end program test_step_cases
