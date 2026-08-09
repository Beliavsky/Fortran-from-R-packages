program roi_example
    use qpoases
    use roi_qpoases
    implicit none

    real(dp) :: q(3,3), linear(3), a(3,3), rhs(3)
    character(len=2) :: direction(3)
    type(qpoases_result) :: result
    integer :: i

    q = 0.0_dp
    do i = 1, 3
        q(i,i) = 1.0_dp
    end do
    linear = [0.0_dp,-5.0_dp,0.0_dp]
    a(1,:) = [-4.0_dp,-3.0_dp,0.0_dp]
    a(2,:) = [ 2.0_dp, 1.0_dp,0.0_dp]
    a(3,:) = [ 0.0_dp,-2.0_dp,1.0_dp]
    rhs = [-8.0_dp,2.0_dp,0.0_dp]
    direction = ">="

    call roi_solve_qp(q,linear,a,direction,rhs,result)

    print '(a,i0)', "status = ", result%status
    print '(a,3f16.10)', "x      = ", result%x
    print '(a,f16.10)', "value  = ", result%objval
end program roi_example
