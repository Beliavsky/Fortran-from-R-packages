program test_assignment
    use lpsolve
    implicit none
    type(lp_result) :: r
    real(dp) :: cost(4,4), expected(4,4), got(4,4)
    integer :: i, j, k

    cost = reshape([2.0_dp,7.0_dp,7.0_dp,2.0_dp, &
                    7.0_dp,7.0_dp,3.0_dp,2.0_dp, &
                    7.0_dp,2.0_dp,8.0_dp,10.0_dp, &
                    1.0_dp,9.0_dp,8.0_dp,2.0_dp], [4,4])
    call lp_assign(cost, r)
    if (r%status /= LP_OPTIMAL) error stop 'assignment status'
    if (abs(r%objective - 8.0_dp) > 1.0e-9_dp) error stop 'assignment objective'

    expected = 0.0_dp
    expected(1,4) = 1.0_dp
    expected(2,3) = 1.0_dp
    expected(3,2) = 1.0_dp
    expected(4,1) = 1.0_dp
    k = 0
    do i = 1, 4
        do j = 1, 4
            k = k + 1
            got(i,j) = r%solution(k)
        end do
    end do
    if (maxval(abs(got - expected)) > 1.0e-9_dp) error stop 'assignment solution'
end program test_assignment
