program basic
    use sadists
    implicit none
    real(dp) :: x(5)

    print '(a,f12.8)', 'doubly noncentral F density: ', &
        ddnf(1.1_dp, 50.0_dp, 80.0_dp, 2.5_dp, 3.5_dp)
    print '(a,f12.8)', 'lambda-prime 0.95 quantile: ', &
        qlambdap(0.95_dp, 50.0_dp, 1.5_dp)

    call rupsilon(x, [30.0_dp, 50.0_dp], [-0.5_dp, 1.0_dp])
    print '(a,5f10.5)', 'five upsilon draws: ', x
end program basic
