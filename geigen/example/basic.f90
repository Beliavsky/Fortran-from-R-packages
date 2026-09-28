program basic
    use geigen_api, only: dp, geigen_result, geigen
    implicit none
    real(dp) :: a(2,2), b(2,2)
    type(geigen_result) :: result

    a = reshape([3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2,2])
    b = 0.0_dp
    b(1,1) = 1.0_dp
    b(2,2) = 2.0_dp
    call geigen(a, b, result)
    if (result%info /= 0) error stop 'geigen failed'
    print '(a,2(1x,f10.6))', 'generalized eigenvalues:', real(result%values, dp)
end program basic
