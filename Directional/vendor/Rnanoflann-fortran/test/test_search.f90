program test_search
    use rnanoflann, only: dp, nn_result, nn
    implicit none
    real(dp) :: data(4,2), points(2,2)
    type(nn_result) :: res

    data = reshape([ &
        0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp, &
        0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], shape(data))
    points = reshape([0.9_dp, 0.0_dp, 0.0_dp, 1.9_dp], shape(points))

    res = nn(data, points, k=2, method="euclidean")
    if (any(shape(res%indices) /= [2,2])) error stop 1
    if (any(res%indices(1,:) /= [2,1])) error stop 2
    if (any(res%indices(2,:) /= [3,1])) error stop 3
    call close(res%distances(1,1), 0.1_dp, 1.0e-14_dp)
    call close(res%distances(2,1), 0.1_dp, 1.0e-14_dp)

    res = nn(data, points, k=2, method="euclidean", square=.true., trans=.false.)
    if (any(shape(res%indices) /= [2,2])) error stop 4
    if (any(res%indices(:,1) /= [2,1])) error stop 5
    call close(res%distances(1,1), 0.01_dp, 1.0e-14_dp)

    print '(a)', 'test_search: PASS'
contains
    subroutine close(value, expected, tol)
        real(dp), intent(in) :: value, expected, tol
        if (abs(value - expected) > tol) error stop 10
    end subroutine close
end program test_search
