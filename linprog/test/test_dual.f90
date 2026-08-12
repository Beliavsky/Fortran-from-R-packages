program test_dual
    use linprog
    implicit none
    real(dp) :: c(3), b(3), a(3,3)
    type(linprog_result) :: r
    type(linprog_control) :: ctl

    c = [1800.0_dp,600.0_dp,600.0_dp]
    b = [40.0_dp,90.0_dp,2500.0_dp]
    a = reshape([0.7_dp,1.5_dp,50.0_dp, 0.35_dp,1.0_dp,12.5_dp, &
                 0.0_dp,3.0_dp,20.0_dp], [3,3])
    ctl%maximum = .true.
    ctl%solve_dual = .true.
    call solveLP(c, b, a, r, ctl)
    call check(r)
    ctl%use_lpsolve = .true.
    call solveLP(c, b, a, r, ctl)
    call check(r)
    print *, 'test_dual: PASS'
contains
    subroutine check(x)
        type(linprog_result), intent(in) :: x
        if (x%status /= 0 .or. x%dual_status /= 0) error stop 'dual status'
        if (abs(x%opt-93600.0_dp) > 1.0e-7_dp) error stop 'dual objective'
        if (maxval(abs(x%con_dual-[0.0_dp,240.0_dp,28.8_dp])) > 1.0e-7_dp) then
            error stop 'dual solution'
        end if
    end subroutine
end program test_dual
