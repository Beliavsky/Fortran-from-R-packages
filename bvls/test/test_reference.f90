program test_reference
    use bvls, only : dp, bvls_fit, bvls_result, bvls_success
    implicit none
    real(dp) :: a(5,3), b(5), lower(3), upper(3)
    real(dp), parameter :: ref(3) = [ &
        2.1232098204137298e-1_dp, 0.0_dp, 7.4380541032052738e-1_dp]
    type(bvls_result) :: fit

    a = reshape([ &
        1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp, &
        2.0_dp,-1.0_dp,0.5_dp,1.5_dp,-2.0_dp, &
        -1.0_dp,0.25_dp,2.0_dp,-0.5_dp,1.0_dp], [5,3])
    b = [1.0_dp,0.0_dp,2.0_dp,-1.0_dp,3.0_dp]
    lower = [-0.2_dp,0.0_dp,-0.5_dp]
    upper = [0.4_dp,1.5_dp,0.8_dp]

    call bvls_fit(a, b, lower, upper, fit)
    if (fit%status /= bvls_success) error stop 'reference: status'
    if (maxval(abs(fit%x-ref)) > 2.0e-13_dp) error stop 'reference: x'
    if (abs(fit%residual_norm-2.5186280305300857_dp) > 2.0e-13_dp) &
        error stop 'reference: norm'
    if (any(fit%istate /= [-2,1,3,1])) error stop 'reference: state'
end program test_reference
