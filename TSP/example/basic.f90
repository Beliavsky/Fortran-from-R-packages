program basic
    use tsp
    implicit none
    real(dp) :: cost(5,5)
    type(tsp_tour) :: tour

    cost = reshape([ &
        0.0_dp, 2.0_dp, 9.0_dp,10.0_dp, 7.0_dp, &
        2.0_dp, 0.0_dp, 6.0_dp, 4.0_dp, 3.0_dp, &
        9.0_dp, 6.0_dp, 0.0_dp, 8.0_dp, 5.0_dp, &
       10.0_dp, 4.0_dp, 8.0_dp, 0.0_dp, 6.0_dp, &
        7.0_dp, 3.0_dp, 5.0_dp, 6.0_dp, 0.0_dp ], [5,5])

    tour = solve_tsp(cost)
    print '(A,*(I0,1X))', 'tour: ', tour%order
    print '(A,F10.3)', 'length: ', tour%length
    print '(A,A)', 'method: ', trim(tour%method)
end program basic
