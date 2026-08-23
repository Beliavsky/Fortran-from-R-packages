program test_parallel
    use rnanoflann, only: dp, nn_result, nn
    implicit none
    real(dp) :: data(5,2), points(3,2)
    type(nn_result) :: a, b

    data(1,:) = [0.0_dp, 0.0_dp]
    data(2,:) = [1.0_dp, 0.0_dp]
    data(3,:) = [0.0_dp, 1.0_dp]
    data(4,:) = [1.0_dp, 1.0_dp]
    data(5,:) = [2.0_dp, 2.0_dp]
    points(1,:) = [0.2_dp, 0.1_dp]
    points(2,:) = [0.8_dp, 0.9_dp]
    points(3,:) = [1.9_dp, 1.8_dp]

    a = nn(data, points, k=3, method="euclidean", parallel=.false.)
    b = nn(data, points, k=3, method="euclidean", parallel=.true., cores=2)
    if (any(a%indices /= b%indices)) error stop 1
    if (maxval(abs(a%distances - b%distances)) > 1.0e-14_dp) error stop 2
    print '(a)', 'test_parallel: PASS'
end program test_parallel
