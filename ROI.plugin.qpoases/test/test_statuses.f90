program test_statuses
    use qpoases
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(dp) :: h(1,1), g(1), a(2,1), lb(1), ub(1), lba(2), uba(2), inf
    type(qpoases_result) :: r

    inf = ieee_value(1.0_dp,ieee_positive_inf)
    h = 1.0_dp
    g = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,1) = 1.0_dp
    lb = -inf
    ub = inf
    lba = [1.0_dp,-inf]
    uba = [inf,0.0_dp]
    call solve_qproblem(h,g,a,lb,ub,lba,uba,r,hessian_type=hst_identity,max_nwsr=500)
    call check(r%status == ret_qp_infeasible, "infeasible status")
    call check(r%infeasible, "infeasible flag")

    h = 0.0_dp
    g = -1.0_dp
    lb = 0.0_dp
    ub = inf
    call solve_qproblemb(h,g,lb,ub,r,hessian_type=hst_zero,max_nwsr=500)
    call check(r%status == ret_qp_unbounded, "unbounded status")
    call check(r%unbounded, "unbounded flag")
    print *, "PASS test_statuses"
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            print *, "FAIL: ", trim(msg), " status=", r%status
            error stop 1
        end if
    end subroutine check
end program test_statuses
