program test_sqproblem_hotstart
    use qpoases
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(qpoases_model) :: model
    real(dp) :: h(2,2), g(2), a(1,2), lb(2), ub(2), lba(1), uba(1)
    real(dp), allocatable :: x(:)
    integer :: status

    h = 0.0_dp
    h(1,1) = 1.0_dp
    h(2,2) = 1.0_dp
    a(1,:) = [1.0_dp,1.0_dp]
    lb = 0.0_dp
    ub = 10.0_dp
    lba = ieee_value(1.0_dp,ieee_negative_inf)
    uba = 1.0_dp
    g = [-2.0_dp,-1.0_dp]
    call init_sqproblem(model,h,g,a,lb,ub,lba,uba,500,status,hessian_type=hst_identity)
    call check(status == ret_qp_solved, "SQ init")
    call get_primal_solution(model,x)
    call check(maxval(abs(x-[1.0_dp,0.0_dp])) < 2.0e-8_dp, "SQ first solution")

    g = [-1.0_dp,-2.0_dp]
    call hotstart_sqproblem(model,h,g,a,lb,ub,lba,uba,500,status)
    call check(status == ret_qp_solved, "SQ hotstart")
    call get_primal_solution(model,x)
    call check(maxval(abs(x-[0.0_dp,1.0_dp])) < 2.0e-8_dp, "SQ hotstart solution")
    print *, "PASS test_sqproblem_hotstart"
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            print *, "FAIL: ", trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_sqproblem_hotstart
