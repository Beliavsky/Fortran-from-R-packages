program test_general_constraints
    use qpoases
    use, intrinsic :: ieee_arithmetic
    implicit none
    real(dp) :: h(2,2), g(2), a(1,2), lb(2), ub(2), lba(1), uba(1), inf
    type(qpoases_result) :: r
    integer :: i

    inf = ieee_value(1.0_dp,ieee_positive_inf)
    h = 0.0_dp
    do i = 1, 2
        h(i,i) = 1.0_dp
    end do
    g = [-1.0_dp,-2.0_dp]
    a(1,:) = [1.0_dp,1.0_dp]
    lb = -inf
    ub = inf
    lba = 1.0_dp
    uba = 1.0_dp
    call solve_qproblem(h,g,a,lb,ub,lba,uba,r,hessian_type=hst_identity)
    call check(r%status == ret_qp_solved, "equality status")
    call check(maxval(abs(r%x-[0.0_dp,1.0_dp])) < 1.0e-8_dp, "equality solution")

    g = [-1.0_dp,-1.0_dp]
    lba = 0.2_dp
    uba = 0.8_dp
    call solve_qproblem(h,g,a,lb,ub,lba,uba,r,hessian_type=hst_identity)
    call check(r%status == ret_qp_solved, "two-sided status")
    call check(maxval(abs(r%x-[0.4_dp,0.4_dp])) < 1.0e-8_dp, "two-sided solution")
    print *, "PASS test_general_constraints"
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            print *, "FAIL: ", trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_general_constraints
