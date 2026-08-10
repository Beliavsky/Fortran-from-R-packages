program test_infinite_bounds
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use bvls, only : dp, bvls_fit, bvls_result, bvls_success
    implicit none
    integer :: i
    real(dp) :: a(20,3), b(20), lower(3), upper(3), truth(3), t, inf
    type(bvls_result) :: fit

    inf = ieee_value(0.0_dp, ieee_positive_inf)
    do i = 1, 20
        t = real(i-1,dp)/10.0_dp
        a(i,1) = exp(-0.5_dp*t)
        a(i,2) = exp(-0.6_dp*t)
        a(i,3) = exp(-t)
    end do
    truth = [1.2_dp,0.5_dp,-0.3_dp]
    b = matmul(a, truth)
    lower = [0.0_dp,0.0_dp,-0.75_dp]
    upper = [inf,inf,0.75_dp]

    call bvls_fit(a, b, lower, upper, fit)
    if (fit%status /= bvls_success) error stop 'infinite: status'
    if (maxval(abs(fit%x-truth)) > 5.0e-12_dp) error stop 'infinite: x'
end program test_infinite_bounds
