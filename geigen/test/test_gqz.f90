program test_gqz
    use geigen_api, only: dp, gqz_result, gqz, gevalues
    implicit none
    type(gqz_result) :: r
    real(dp) :: a(3,3), b(3,3), err
    complex(dp) :: ac(2,2), bc(2,2)
    complex(dp), allocatable :: vals(:)

    a = reshape([1.0_dp, 4.0_dp, 0.0_dp, &
                 2.0_dp, 3.0_dp, 1.0_dp, &
                 0.0_dp, 1.0_dp, -2.0_dp], [3,3])
    b = 0.0_dp
    b(1,1) = 1.0_dp
    b(2,2) = 2.0_dp
    b(3,3) = 1.5_dp
    call gqz(a, b, r)
    if (r%info /= 0) error stop 'real gqz status'
    err = maxval(abs(cmplx(a,0.0_dp,dp) - matmul(r%q, matmul(r%s, transpose(r%z)))))
    if (err > 1.0e-10_dp) error stop 'real gqz A reconstruction'
    err = maxval(abs(cmplx(b,0.0_dp,dp) - matmul(r%q, matmul(r%t, transpose(r%z)))))
    if (err > 1.0e-10_dp) error stop 'real gqz B reconstruction'
    vals = gevalues(r)
    if (maxval(abs(vals-r%values)) > 1.0e-14_dp) error stop 'gevalues real'

    call gqz(a, b, r, '-')
    if (r%info /= 0) error stop 'real ordered gqz status'
    if (r%sdim < 1) error stop 'real ordered gqz count'

    ac = reshape([cmplx(1.0_dp,1.0_dp,dp), cmplx(2.0_dp,-1.0_dp,dp), &
                  cmplx(0.5_dp,0.2_dp,dp), cmplx(-1.0_dp,0.5_dp,dp)], [2,2])
    bc = cmplx(0.0_dp,0.0_dp,dp)
    bc(1,1) = cmplx(1.0_dp,0.0_dp,dp)
    bc(2,2) = cmplx(2.0_dp,0.0_dp,dp)
    call gqz(ac, bc, r)
    if (r%info /= 0) error stop 'complex gqz status'
    err = maxval(abs(ac - matmul(r%q, matmul(r%s, conjg(transpose(r%z))))))
    if (err > 1.0e-10_dp) error stop 'complex gqz A reconstruction'
    err = maxval(abs(bc - matmul(r%q, matmul(r%t, conjg(transpose(r%z))))))
    if (err > 1.0e-10_dp) error stop 'complex gqz B reconstruction'

    print *, 'test_gqz: PASS'
end program test_gqz
