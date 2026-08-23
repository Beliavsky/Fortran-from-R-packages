program basic
    use rnanoflann, only: dp, nn_result, nn
    implicit none
    real(dp) :: data(5,2), query(2,2)
    type(nn_result) :: ans
    integer :: i

    data(1,:) = [0.0_dp, 0.0_dp]
    data(2,:) = [1.0_dp, 0.0_dp]
    data(3,:) = [0.0_dp, 1.0_dp]
    data(4,:) = [1.0_dp, 1.0_dp]
    data(5,:) = [2.0_dp, 2.0_dp]
    query(1,:) = [0.9_dp, 0.1_dp]
    query(2,:) = [1.8_dp, 1.9_dp]

    ans = nn(data, query, k=2, method="euclidean")
    do i = 1, size(query,1)
        print '(a,i0,a,2(i0,1x))', 'query ', i, ' neighbors: ', ans%indices(i,:)
        print '(a,2(f10.6,1x))', 'distances: ', ans%distances(i,:)
    end do
end program basic
