program test_warm_start
    use bvls, only : dp, bvls_fit, bvls_result, bvls_success
    implicit none
    real(dp) :: a(6,3), b(6), lower(3), upper(3)
    type(bvls_result) :: cold, warm

    a = reshape([ &
        1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
        1.0_dp,-1.0_dp,2.0_dp,-2.0_dp,0.5_dp,1.5_dp, &
        0.5_dp,1.0_dp,-0.5_dp,2.0_dp,-1.0_dp,0.25_dp], [6,3])
    b = [1.0_dp,0.5_dp,2.0_dp,-1.0_dp,2.5_dp,0.0_dp]
    lower = [-0.5_dp,0.0_dp,-0.25_dp]
    upper = [0.75_dp,2.0_dp,1.25_dp]

    call bvls_fit(a, b, lower, upper, cold)
    call bvls_fit(a, b, lower, upper, warm, key=1, istate=cold%istate)
    if (cold%status /= bvls_success .or. warm%status /= bvls_success) &
        error stop 'warm: status'
    if (maxval(abs(cold%x-warm%x)) > 5.0e-12_dp) error stop 'warm: x'
    if (abs(cold%deviance-warm%deviance) > 5.0e-12_dp) error stop 'warm: rss'
end program test_warm_start
