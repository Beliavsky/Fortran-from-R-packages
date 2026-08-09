program test_qproblemb_hotstart
    use qpoases
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(qpoases_model) :: model
    real(dp) :: h(2,2), g(2), lb(2), ub(2)
    real(dp), allocatable :: x(:)
    integer :: status

    h = 0.0_dp
    h(1,1) = 1.0_dp
    h(2,2) = 1.0_dp
    g = [-1.0_dp,-2.0_dp]
    lb = 0.0_dp
    ub = 10.0_dp
    call init_qproblemb(model,h,g,lb,ub,200,status,hessian_type=hst_identity)
    call check(status == ret_qp_solved, "initial solve")
    call get_primal_solution(model,x)
    call check(maxval(abs(x-[1.0_dp,2.0_dp])) < 1.0e-9_dp, "initial solution")
    call check(is_initialised(model) .and. is_solved(model), "model state")
    call check(get_number_of_variables(model) == 2, "variable count")

    g = [-3.0_dp,-4.0_dp]
    call hotstart_qproblemb(model,g,lb,ub,200,status)
    call check(status == ret_qp_solved, "hotstart status")
    call get_primal_solution(model,x)
    call check(maxval(abs(x-[3.0_dp,4.0_dp])) < 1.0e-8_dp, "hotstart solution")
    print *, "PASS test_qproblemb_hotstart"
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            print *, "FAIL: ", trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_qproblemb_hotstart
