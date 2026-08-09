program test_roi_regressions
    use qpoases
    use roi_qpoases
    implicit none
    real(dp) :: q3(3,3), l3(3), a3(3,3), rhs3(3)
    real(dp) :: q2(2,2), l2(2), a2(2,2), rhs2(2), tiny
    character(len=2) :: d3(3), d2(2)
    type(qpoases_result) :: r
    integer :: i

    q3 = 0.0_dp
    do i = 1, 3
        q3(i,i) = 1.0_dp
    end do
    l3 = [0.0_dp,-5.0_dp,0.0_dp]
    a3(1,:) = [-4.0_dp,-3.0_dp,0.0_dp]
    a3(2,:) = [ 2.0_dp, 1.0_dp,0.0_dp]
    a3(3,:) = [ 0.0_dp,-2.0_dp,1.0_dp]
    rhs3 = [-8.0_dp,2.0_dp,0.0_dp]
    d3 = ">="
    call roi_solve_qp(q3,l3,a3,d3,rhs3,r)
    call check(r%status == ret_qp_solved, "ROI QP 1 status")
    call check(maxval(abs(r%x-[0.476190476190476_dp,1.04761904761905_dp, &
        2.09523809523810_dp])) < 2.0e-8_dp, "ROI QP 1 solution")
    call check(abs(r%objval + 2.38095238095238_dp) < 2.0e-8_dp, "ROI QP 1 objective")

    tiny = epsilon(1.0_dp) * 100.0_dp
    q2 = 0.0_dp
    q2(1,1) = 1.0_dp
    q2(2,2) = tiny
    l2 = [-2.0_dp,1.0_dp]
    a2(1,:) = [1.0_dp,0.0_dp]
    a2(2,:) = [1.0_dp,0.0_dp]
    rhs2 = [3.0_dp,0.0_dp]
    d2 = ["<=",">="]
    call roi_solve_qp(q2,l2,a2,d2,rhs2,r)
    call check(maxval(abs(r%x-[2.0_dp,0.0_dp])) < 1.0e-8_dp, "ROI QP 2 solution")
    call check(abs(r%objval + 2.0_dp) < 1.0e-8_dp, "ROI QP 2 objective")

    call roi_solve_qp(-q3,-l3,a3,d3,rhs3,r,maximum=.true.)
    call check(maxval(abs(r%x-[0.476190476190476_dp,1.04761904761905_dp, &
        2.09523809523810_dp])) < 2.0e-8_dp, "ROI maximize solution")
    call check(abs(r%objval - 2.38095238095238_dp) < 2.0e-8_dp, "ROI maximize objective")
    print *, "PASS test_roi_regressions"
contains
    subroutine check(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            print *, "FAIL: ", trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_roi_regressions
