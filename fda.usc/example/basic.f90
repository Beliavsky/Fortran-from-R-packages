program basic_fdausc
    use fdausc
    implicit none

    real(dp) :: t(5)
    real(dp) :: curves(3, 5)
    real(dp) :: distances(3, 3)

    t = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
    curves(1, :) = t
    curves(2, :) = 2.0_dp*t
    curves(3, :) = 3.0_dp*t
    call lp_distance_matrix(t, curves, curves, distances, p=2.0_dp, method=2)
    print '(a,f10.6)', "L2 distance between curves 1 and 3: ", distances(1, 3)
end program basic_fdausc
