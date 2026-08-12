program assignment_example
    use lpsolve
    implicit none
    type(lp_result) :: result
    real(dp) :: cost(4,4), assignment(4,4)
    integer :: i

    cost = reshape([2.0_dp,7.0_dp,7.0_dp,2.0_dp, &
                    7.0_dp,7.0_dp,3.0_dp,2.0_dp, &
                    7.0_dp,2.0_dp,8.0_dp,10.0_dp, &
                    1.0_dp,9.0_dp,8.0_dp,2.0_dp], [4,4])

    call lp_assign(cost, result, assignment=assignment)
    print '(a,f10.3)', 'objective: ', result%objective
    do i = 1, 4
        print '(4f6.1)', assignment(i,:)
    end do
end program assignment_example
