program test_lindsey_optimizer
    use multiplicative_multinomial
    implicit none
    type(paras_type) :: par
    type(glm_fit_type) :: gl
    type(mm_fit_type) :: fit
    real(dp), allocatable :: prob(:), th(:,:)
    real(dp) :: wt(3)
    integer :: obs(3,2)

    obs(1,:) = [2, 0]
    obs(2,:) = [1, 1]
    obs(3,:) = [0, 2]
    wt = [16.0_dp, 48.0_dp, 36.0_dp]
    call lindsey_fit(obs, par, gl, wt)
    if (.not. gl%converged) error stop "Lindsey Poisson GLM did not converge"
    prob = p(par)
    th = theta(par)
    call check_close(prob(1), 0.4_dp, 2.0e-9_dp, "Lindsey p1")
    call check_close(prob(2), 0.6_dp, 2.0e-9_dp, "Lindsey p2")
    call check_close(th(1,2), 1.0_dp, 2.0e-8_dp, "Lindsey theta")

    call optimizer_allsamesum(obs, fit, wt, par, max_iter=500, tol=1.0e-10_dp)
    prob = p(fit%parameters)
    th = theta(fit%parameters)
    call check_close(prob(1), 0.4_dp, 3.0e-6_dp, "optimizer p1")
    call check_close(th(1,2), 1.0_dp, 2.0e-5_dp, "optimizer theta")
    call check_close(fit%loglik, mm_loglik(obs, par, wt), 2.0e-8_dp, "optimizer loglik")

    call optimizer_allsamesum(obs, fit, wt, par, max_iter=500, tol=1.0e-10_dp, method="Nelder")
    prob = p(fit%parameters)
    th = theta(fit%parameters)
    call check_close(prob(1), 0.4_dp, 3.0e-6_dp, "Nelder p1")
    call check_close(th(1,2), 1.0_dp, 2.0e-5_dp, "Nelder theta")

    print '(a)', 'test_lindsey_optimizer: PASS'

contains

    subroutine check_close(x, ref, tol, label)
        real(dp), intent(in) :: x, ref, tol
        character(len=*), intent(in) :: label
        if (abs(x-ref) > tol * max(1.0_dp, abs(ref))) then
            write(*,'(a,2es24.14)') trim(label)//' failed: ', x, ref
            error stop 1
        end if
    end subroutine check_close

end program test_lindsey_optimizer
