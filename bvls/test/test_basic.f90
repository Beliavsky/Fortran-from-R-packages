program test_basic
    use bvls, only : dp, bvls_fit, bvls_result, bvls_success
    implicit none
    real(dp) :: a(3,2), b(3), lower(2), upper(2)
    type(bvls_result) :: fit

    a = reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,1.0_dp], [3,2])
    b = [2.0_dp,-1.0_dp,0.5_dp]
    lower = [0.0_dp,0.0_dp]
    upper = [1.0_dp,2.0_dp]

    call bvls_fit(a, b, lower, upper, fit)
    if (fit%status /= bvls_success) error stop 'basic: status'
    if (maxval(abs(fit%x-[1.0_dp,0.0_dp])) > 1.0e-12_dp) error stop 'basic: x'
    if (abs(fit%deviance-2.25_dp) > 1.0e-12_dp) error stop 'basic: deviance'
    if (abs(fit%residual_norm-1.5_dp) > 1.0e-12_dp) error stop 'basic: norm'
end program test_basic
