program test_radius
    use rnanoflann, only: dp, nn_result, nn
    implicit none
    real(dp) :: data(4,2), points(2,2), missing
    type(nn_result) :: res

    data = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp, &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], shape(data))
    points = reshape([0.9_dp, 0.0_dp, 0.0_dp, 1.9_dp], shape(points))
    missing = sqrt(huge(0.0_dp))

    res = nn(data, points, k=2, method="euclidean", search="radius", radius=0.25_dp, sorted=.true.)
    if (any(res%counts /= [1,1])) error stop 1
    if (res%indices(1,1) /= 2 .or. res%indices(2,1) /= 3) error stop 2
    if (res%indices(1,2) /= 0 .or. res%indices(2,2) /= 0) error stop 3
    if (res%distances(1,2) < 0.99_dp * missing .or. res%distances(2,2) < 0.99_dp * missing) error stop 4

    res = nn(data, points(1:1,:), k=2, method="euclidean", search="radius", radius=2.20_dp, sorted=.true.)
    if (res%counts(1) /= 4) error stop 5
    if (any(res%indices(1,:) /= [2,1])) error stop 6

    print '(a)', 'test_radius: PASS'
end program test_radius
