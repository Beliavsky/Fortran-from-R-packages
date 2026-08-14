program basic_network
    use sna, only : dp, geodist_result, geodist, degree, betweenness, &
                    dyad_census, triad_census, gden
    implicit none

    real(dp) :: g(4,4), dc(3)
    real(dp), allocatable :: deg(:), bet(:), tc(:)
    type(geodist_result) :: gd

    g = 0.0_dp
    g(1,2) = 1.0_dp
    g(2,3) = 1.0_dp
    g(3,1) = 1.0_dp
    g(3,4) = 1.0_dp

    gd = geodist(g)
    deg = degree(g, cmode='outdegree', ignore_eval=.true.)
    bet = betweenness(g)
    dc = dyad_census(g)
    tc = triad_census(g)

    print '(a,f8.4)', 'density: ', gden(g)
    print '(a,*(f8.3,1x))', 'outdegree: ', deg
    print '(a,*(f8.3,1x))', 'betweenness: ', bet
    print '(a,*(f8.1,1x))', 'dyad census (M,A,N): ', dc
    print '(a,f8.1)', 'number of triads: ', sum(tc)
    print '(a,f8.1)', 'distance 1 -> 4: ', gd%distance(1,4)
end program basic_network
