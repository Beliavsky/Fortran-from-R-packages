program test_strategies
    use deoptim, only : dp, de_control, de_result, deoptim_solve, de_success
    implicit none
    type(de_control) :: ctrl
    type(de_result) :: res
    real(dp) :: lo(2), hi(2)
    integer :: strategy

    lo = -5.0_dp
    hi = 5.0_dp

    do strategy = 1, 6
        ctrl = de_control()
        ctrl%strategy = strategy
        ctrl%np = 60
        ctrl%itermax = 700
        ctrl%f = 0.8_dp
        ctrl%cr = 0.9_dp
        ctrl%seed = 1000 + strategy
        ctrl%reltol = 0.0_dp
        ctrl%steptol = ctrl%itermax
        if (strategy == 6) then
            ctrl%p = 0.2_dp
            ctrl%c = 0.1_dp
        end if
        call deoptim_solve(rosenbrock, lo, hi, res, ctrl)
        if (res%status /= de_success) error stop "strategy solve failed"
        if (res%bestval > 1.0e-7_dp) then
            write(*,*) "strategy", strategy, "best", res%bestval, res%bestmem
            error stop "Rosenbrock convergence failure"
        end if
        if (res%nfeval /= int(ctrl%np * (res%iter + 1), kind(res%nfeval))) &
            error stop "evaluation count mismatch"
    end do

    print *, "test_strategies: PASS"
contains
    function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 100.0_dp * (x(2) - x(1) * x(1))**2 + (1.0_dp - x(1))**2
    end function rosenbrock
end program test_strategies
