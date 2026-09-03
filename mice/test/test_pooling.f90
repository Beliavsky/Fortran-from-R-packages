program test_pooling
    use mice, only : dp, mice_ok, pool_scalar_result, pool_vector_result, pool_scalar, pool_vector, pooled_wald
    implicit none
    real(dp) :: q(2), u(2), qv(2, 3), uv(2, 2, 3), stat
    type(pool_scalar_result) :: scalar
    type(pool_vector_result) :: vector
    integer :: info, i

    q = [-1.5457_dp, -1.428_dp]
    u = [0.9723_dp**2, 1.041_dp**2]
    call pool_scalar(q, u, scalar, info, n=25.0_dp, k=2)
    if (info /= mice_ok) error stop "pool scalar status"
    if (abs(scalar%qbar + 1.48685_dp) > 1.0e-10_dp) error stop "qbar reference"
    if (scalar%total <= scalar%ubar) error stop "Rubin total"
    qv(:, 1) = [1.0_dp, 2.0_dp]
    qv(:, 2) = [1.2_dp, 1.8_dp]
    qv(:, 3) = [0.8_dp, 2.2_dp]
    uv = 0.0_dp
    do i = 1, 3
        uv(1, 1, i) = 0.04_dp
        uv(2, 2, i) = 0.09_dp
    end do
    call pool_vector(qv, uv, vector, info)
    if (info /= mice_ok) error stop "pool vector status"
    if (maxval(abs(vector%qbar - [1.0_dp, 2.0_dp])) > 1.0e-13_dp) error stop "vector qbar"
    call pooled_wald(vector%qbar, [0.0_dp, 0.0_dp], vector%total, stat, info)
    if (info /= mice_ok .or. stat <= 0.0_dp) error stop "pooled Wald"
    print *, "test_pooling: PASS"
end program test_pooling
