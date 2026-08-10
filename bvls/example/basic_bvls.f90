program basic_bvls
    use bvls, only : dp, bvls_fit, bvls_result
    implicit none
    real(dp) :: a(3,2), b(3), lower(2), upper(2)
    type(bvls_result) :: fit

    a = reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,1.0_dp], [3,2])
    b = [2.0_dp,-1.0_dp,0.5_dp]
    lower = [0.0_dp,0.0_dp]
    upper = [1.0_dp,2.0_dp]

    call bvls_fit(a, b, lower, upper, fit)
    print '(a,*(f12.6,1x))', 'x = ', fit%x
    print '(a,es14.6)', 'RSS = ', fit%deviance
end program basic_bvls
