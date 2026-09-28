program test_gsvd_rank
    use geigen_api, only: dp, gsvd_result, gsvd, gsvd_or, gsvd_d1, gsvd_d2
    implicit none
    type(gsvd_result) :: z
    real(dp) :: a(3,4), b(2,4), err
    complex(dp), allocatable :: ormat(:,:)
    real(dp), allocatable :: d1(:,:), d2(:,:)

    a = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,1) = 1.0_dp
    a(3,2) = -1.0_dp
    b = 0.0_dp
    b(1,1) = 0.5_dp
    b(1,2) = 1.0_dp
    b(2,1) = -1.0_dp
    b(2,2) = 0.25_dp

    call gsvd(a, b, z)
    if (z%info /= 0) error stop 'rank-deficient gsvd status'
    if (z%rank /= 2) error stop 'rank-deficient gsvd rank'
    ormat = gsvd_or(z)
    d1 = gsvd_d1(z)
    d2 = gsvd_d2(z)
    err = maxval(abs(matmul(conjg(transpose(z%u)), matmul(cmplx(a,0.0_dp,dp), z%q)) - &
        matmul(cmplx(d1,0.0_dp,dp), ormat)))
    if (err > 1.0e-9_dp) error stop 'rank-deficient gsvd A identity'
    err = maxval(abs(matmul(conjg(transpose(z%v)), matmul(cmplx(b,0.0_dp,dp), z%q)) - &
        matmul(cmplx(d2,0.0_dp,dp), ormat)))
    if (err > 1.0e-9_dp) error stop 'rank-deficient gsvd B identity'

    print *, 'test_gsvd_rank: PASS'
end program test_gsvd_rank
