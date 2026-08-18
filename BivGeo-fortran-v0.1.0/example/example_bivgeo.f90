program example_bivgeo
    use bivgeo, only : dp, bivgeo_params, make_bivgeo_params, dbivgeo2, pbivgeo, covbivgeo, corbivgeo, cfbivgeo
    implicit none

    type(bivgeo_params) :: theta

    theta = make_bivgeo_params(0.2_dp, 0.4_dp, 0.7_dp)
    print '(a,f14.8)', 'P(X=1,Y=2) = ', dbivgeo2(1, 2, theta)
    print '(a,f14.8)', 'P(X<=1,Y<=2) = ', pbivgeo(1, 2, theta)
    print '(a,f14.8)', 'Cov(X,Y) = ', covbivgeo(theta)
    print '(a,f14.8)', 'Cor(X,Y) = ', corbivgeo(theta)
    print '(a,f14.8)', 'E[X Y] = ', cfbivgeo(theta)
end program example_bivgeo
