program test_geigen
    use geigen_api, only: dp, geigen_result, geigen
    implicit none
    type(geigen_result) :: r
    real(dp) :: a(3,3), b(3,3), err
    complex(dp) :: ac(2,2), bc(2,2)
    integer :: j

    a = reshape([4.0_dp, 1.0_dp, 0.0_dp, &
                 1.0_dp, 3.0_dp, 0.5_dp, &
                 0.0_dp, 0.5_dp, 2.0_dp], [3,3])
    b = 0.0_dp
    b(1,1) = 2.0_dp
    b(2,2) = 1.5_dp
    b(3,3) = 1.0_dp
    call geigen(a, b, r)
    if (r%info /= 0) error stop 'symmetric geigen status'
    if (.not. r%symmetric_used) error stop 'symmetric path not used'
    do j = 1, 3
        err = maxval(abs(matmul(cmplx(a,0.0_dp,dp), r%vectors(:,j)) - &
            r%values(j) * matmul(cmplx(b,0.0_dp,dp), r%vectors(:,j))))
        if (err > 1.0e-10_dp) error stop 'symmetric geigen residual'
    end do

    a = 0.0_dp
    a(1,2) = -1.0_dp
    a(2,1) = 1.0_dp
    a(3,3) = 2.0_dp
    b = 0.0_dp
    b(1,1) = 1.0_dp
    b(2,2) = 1.0_dp
    b(3,3) = 1.0_dp
    call geigen(a, b, r, symmetric=.false.)
    if (r%info /= 0) error stop 'general real geigen status'
    do j = 1, 3
        err = maxval(abs(matmul(cmplx(a,0.0_dp,dp), r%vectors(:,j)) - &
            r%values(j) * matmul(cmplx(b,0.0_dp,dp), r%vectors(:,j))))
        if (err > 1.0e-10_dp) error stop 'general real geigen residual'
    end do

    ac = reshape([cmplx(2.0_dp,0.0_dp,dp), cmplx(1.0_dp,-1.0_dp,dp), &
                  cmplx(1.0_dp,1.0_dp,dp), cmplx(3.0_dp,0.0_dp,dp)], [2,2])
    bc = cmplx(0.0_dp,0.0_dp,dp)
    bc(1,1) = cmplx(2.0_dp,0.0_dp,dp)
    bc(2,2) = cmplx(1.0_dp,0.0_dp,dp)
    call geigen(ac, bc, r)
    if (r%info /= 0) error stop 'Hermitian geigen status'
    do j = 1, 2
        err = maxval(abs(matmul(ac, r%vectors(:,j)) - r%values(j) * matmul(bc, r%vectors(:,j))))
        if (err > 1.0e-10_dp) error stop 'Hermitian geigen residual'
    end do

    print *, 'test_geigen: PASS'
end program test_geigen
