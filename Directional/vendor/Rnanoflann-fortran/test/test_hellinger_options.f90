program test_hellinger_options
    use rnanoflann, only: dp, nn_result, nn
    implicit none
    real(dp) :: data(3,3), points(1,3)
    type(nn_result) :: res

    data(1,:) = [0.1_dp, 0.2_dp, 0.7_dp]
    data(2,:) = [0.2_dp, 0.5_dp, 0.3_dp]
    data(3,:) = [0.4_dp, 0.4_dp, 0.2_dp]
    points(1,:) = [0.2_dp, 0.5_dp, 0.3_dp]

    res = nn(data, points, k=2, method="hellinger")
    if (res%indices(1,1) /= 2) error stop 1
    if (abs(res%distances(1,1)) > 1.0e-14_dp) error stop 2

    ! Upstream calls these similarities "distances" and minimizes them.
    data = 0.0_dp
    data(1,:) = [1.0_dp, 0.0_dp, 0.0_dp]
    data(2,:) = [0.0_dp, 1.0_dp, 0.0_dp]
    data(3,:) = [-1.0_dp, 0.0_dp, 0.0_dp]
    points(1,:) = [1.0_dp, 0.0_dp, 0.0_dp]
    res = nn(data, points, k=1, method="cosine")
    if (res%indices(1,1) /= 3) error stop 3

    print '(a)', 'test_hellinger_options: PASS'
end program test_hellinger_options
