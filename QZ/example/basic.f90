program basic
    use qz_api, only: dp, qz_geigen, qz_eigen_real_result
    implicit none

    real(dp) :: a(2, 2)
    type(qz_eigen_real_result) :: result

    a = reshape([0.0_dp, 1.0_dp, -2.0_dp, -3.0_dp], [2, 2])
    call qz_geigen(a, result, only_values=.true.)
    if (result%info /= 0) error stop 'qz_geigen failed'
    print '(a,2(1x,"(",f8.4,",",f8.4,")"))', 'eigenvalues:', result%w
end program basic
