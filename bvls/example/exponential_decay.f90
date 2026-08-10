program exponential_decay
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use bvls, only : dp, bvls_fit, bvls_result
    implicit none
    integer :: i
    real(dp) :: a(51,3), b(51), lower(3), upper(3), truth(3), t, inf
    type(bvls_result) :: fit

    inf = ieee_value(0.0_dp, ieee_positive_inf)
    do i = 1, 51
        t = 0.04_dp * real(i-1,dp)
        a(i,1) = exp(-0.5_dp*t)
        a(i,2) = exp(-0.6_dp*t)
        a(i,3) = exp(-t)
    end do
    truth = [0.7_dp,0.2_dp,-0.4_dp]
    b = matmul(a, truth)
    lower = [0.0_dp,0.0_dp,-0.75_dp]
    upper = [inf,inf,0.75_dp]

    call bvls_fit(a, b, lower, upper, fit)
    print '(a,*(f12.7,1x))', 'truth = ', truth
    print '(a,*(f12.7,1x))', 'fit   = ', fit%x
    print '(a,es14.6)', 'residual norm = ', fit%residual_norm
end program exponential_decay
