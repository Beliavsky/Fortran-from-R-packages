program test_gsvd
    use geigen_api, only: dp, gsvd_result, gsvd, gsvd_or, gsvd_d1, gsvd_d2, gsvd_r
    implicit none
    type(gsvd_result) :: z
    real(dp) :: a(4,3), b(2,3), ar(2,4), br(3,4), err
    complex(dp) :: ac(3,3), bc(2,3)
    complex(dp), allocatable :: ormat(:,:), rr(:,:)
    real(dp), allocatable :: d1(:,:), d2(:,:)

    a = reshape([1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp, &
                 0.0_dp, 1.0_dp, 2.0_dp, -0.5_dp, &
                 2.0_dp, 0.0_dp, 1.0_dp, 1.5_dp], [4,3])
    b = reshape([1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, -1.0_dp, 2.0_dp], [2,3])
    call gsvd(a, b, z)
    if (z%info /= 0) error stop 'real gsvd status'
    if (z%rank < 1) error stop 'real gsvd rank'
    ormat = gsvd_or(z)
    rr = gsvd_r(z)
    d1 = gsvd_d1(z)
    d2 = gsvd_d2(z)
    if (maxval(abs(rr-z%r)) > 1.0e-14_dp) error stop 'gsvd.R helper'
    err = maxval(abs(matmul(conjg(transpose(z%u)), matmul(cmplx(a,0.0_dp,dp), z%q)) - &
        matmul(cmplx(d1,0.0_dp,dp), ormat)))
    if (err > 5.0e-10_dp) error stop 'real gsvd A identity'
    err = maxval(abs(matmul(conjg(transpose(z%v)), matmul(cmplx(b,0.0_dp,dp), z%q)) - &
        matmul(cmplx(d2,0.0_dp,dp), ormat)))
    if (err > 5.0e-10_dp) error stop 'real gsvd B identity'
    if (maxval(abs(z%alpha(1:z%rank)**2 + z%beta(1:z%rank)**2 - 1.0_dp)) > 1.0e-10_dp) &
        error stop 'real gsvd alpha beta'

    ! Wide pair exercises the leading zero block in [0 R].
    ar = reshape([1.0_dp,0.0_dp, 0.0_dp,1.0_dp, 1.0_dp,1.0_dp, 2.0_dp,-1.0_dp], [2,4])
    br = reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,1.0_dp, 1.0_dp,1.0_dp,0.0_dp, 0.5_dp,-1.0_dp,2.0_dp], [3,4])
    call gsvd(ar, br, z)
    if (z%info /= 0) error stop 'wide real gsvd status'
    ormat = gsvd_or(z)
    d1 = gsvd_d1(z)
    d2 = gsvd_d2(z)
    err = maxval(abs(matmul(conjg(transpose(z%u)), matmul(cmplx(ar,0.0_dp,dp), z%q)) - &
        matmul(cmplx(d1,0.0_dp,dp), ormat)))
    if (err > 1.0e-9_dp) error stop 'wide real gsvd A identity'
    err = maxval(abs(matmul(conjg(transpose(z%v)), matmul(cmplx(br,0.0_dp,dp), z%q)) - &
        matmul(cmplx(d2,0.0_dp,dp), ormat)))
    if (err > 1.0e-9_dp) error stop 'wide real gsvd B identity'

    ac = reshape([cmplx(1.0_dp,0.2_dp,dp), cmplx(0.0_dp,1.0_dp,dp), cmplx(2.0_dp,-0.5_dp,dp), &
                  cmplx(0.5_dp,-0.1_dp,dp), cmplx(1.0_dp,0.0_dp,dp), cmplx(-1.0_dp,0.3_dp,dp), &
                  cmplx(2.0_dp,0.0_dp,dp), cmplx(0.5_dp,0.5_dp,dp), cmplx(1.0_dp,-1.0_dp,dp)], [3,3])
    bc = reshape([cmplx(1.0_dp,0.0_dp,dp), cmplx(0.5_dp,0.2_dp,dp), &
                  cmplx(0.0_dp,1.0_dp,dp), cmplx(1.0_dp,0.0_dp,dp), &
                  cmplx(1.0_dp,-0.2_dp,dp), cmplx(-0.5_dp,0.4_dp,dp)], [2,3])
    call gsvd(ac, bc, z)
    if (z%info /= 0) error stop 'complex gsvd status'
    ormat = gsvd_or(z)
    d1 = gsvd_d1(z)
    d2 = gsvd_d2(z)
    err = maxval(abs(matmul(conjg(transpose(z%u)), matmul(ac, z%q)) - matmul(cmplx(d1,0.0_dp,dp), ormat)))
    if (err > 1.0e-9_dp) error stop 'complex gsvd A identity'
    err = maxval(abs(matmul(conjg(transpose(z%v)), matmul(bc, z%q)) - matmul(cmplx(d2,0.0_dp,dp), ormat)))
    if (err > 1.0e-9_dp) error stop 'complex gsvd B identity'

    print *, 'test_gsvd: PASS'
end program test_gsvd
