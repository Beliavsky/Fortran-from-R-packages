program test_model_diagnostics
    use qpoases
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(qpoases_model) :: model
    real(dp) :: h(2,2), g(2), a(1,2), lb(2), ub(2), lba(1), uba(1), inf
    integer :: status

    inf = ieee_value(1.0_dp,ieee_positive_inf)
    h = 0.0_dp
    h(1,1) = 1.0_dp
    h(2,2) = 1.0_dp
    g = [-2.0_dp,-2.0_dp]
    a(1,:) = [1.0_dp,1.0_dp]
    lb = 0.0_dp
    ub = inf
    lba = -inf
    uba = 1.0_dp
    call init_qproblem(model,h,g,a,lb,ub,lba,uba,500,status,hessian_type=hst_identity)
    call check(status == ret_qp_solved, "diagnostic solve")
    call check(get_number_of_constraints(model) == 1, "constraint count")
    call check(get_number_of_active_constraints(model) == 1, "active constraint count")
    call check(get_number_of_inactive_constraints(model) == 0, "inactive constraint count")
    call check(get_number_of_equality_constraints(model) == 0, "equality count")
    call check(abs(get_objval(model) + 1.75_dp) < 1.0e-8_dp, "objective getter")
    print *, "PASS test_model_diagnostics"
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            print *, "FAIL: ", trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_model_diagnostics
