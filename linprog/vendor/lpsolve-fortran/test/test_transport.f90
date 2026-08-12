program test_transport
    use lpsolve
    implicit none
    type(lp_result) :: r
    real(dp) :: cost(8,5), row_rhs(8), col_rhs(5), got(8,5), expected(8,5)
    integer :: row_sense(8), col_sense(5), i, j, k

    cost = 10000.0_dp
    cost(4,1) = 0.0_dp
    cost(1:3,5) = 0.0_dp
    cost(5:8,5) = 0.0_dp
    cost(1,2) = 7.0_dp
    cost(2,3) = 7.0_dp
    cost(3,4) = 7.0_dp
    cost(1,3) = 7.7_dp
    cost(2,4) = 7.7_dp
    cost(5,1) = 8.0_dp
    cost(7,3) = 8.0_dp
    cost(1,4) = 8.4_dp
    cost(6,2) = 9.0_dp
    cost(8,4) = 10.0_dp
    cost(4,2:4) = [0.7_dp,1.4_dp,2.1_dp]

    row_rhs = [200.0_dp,300.0_dp,350.0_dp,200.0_dp,100.0_dp,50.0_dp,100.0_dp,150.0_dp]
    col_rhs = [250.0_dp,100.0_dp,400.0_dp,500.0_dp,200.0_dp]
    row_sense = LP_LE
    col_sense = LP_GE
    call lp_transport(cost, row_sense, row_rhs, col_sense, col_rhs, r)
    if (r%status /= LP_OPTIMAL) error stop 'transport status'
    if (abs(r%objective - 7790.0_dp) > 1.0e-6_dp) error stop 'transport objective'

    expected = 0.0_dp
    expected(1,2) = 100.0_dp
    expected(1,4) = 100.0_dp
    expected(2,3) = 300.0_dp
    expected(3,4) = 350.0_dp
    expected(4,1) = 200.0_dp
    expected(5,1) = 50.0_dp
    expected(5,5) = 50.0_dp
    expected(6,5) = 50.0_dp
    expected(7,3) = 100.0_dp
    expected(8,4) = 50.0_dp
    expected(8,5) = 100.0_dp
    k = 0
    do i = 1, 8
        do j = 1, 5
            k = k + 1
            got(i,j) = r%solution(k)
        end do
    end do
    if (maxval(abs(got - expected)) > 1.0e-6_dp) error stop 'transport solution'
end program test_transport
