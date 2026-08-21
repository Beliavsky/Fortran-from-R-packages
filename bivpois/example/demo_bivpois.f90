program demo_bivpois
    use bivpois, only : dp, seed_rng, rbp, bp_mle2_result, bp_mle2, dbp_scalar
    implicit none
    integer, parameter :: n = 300
    integer :: x(n, 2)
    real(dp), parameter :: lambda_true(3) = [3.0_dp, 5.0_dp, 2.0_dp]
    type(bp_mle2_result) :: fit

    call seed_rng(12345)
    call rbp(n, lambda_true, x)
    fit = bp_mle2(x(:, 1), x(:, 2))

    print '(a,3f12.6)', 'true lambda:      ', lambda_true
    print '(a,3f12.6)', 'estimated lambda: ', fit%lambda
    print '(a,f14.8)', 'P(X=3,Y=4):     ', dbp_scalar(3, 4, lambda_true, .false.)
end program demo_bivpois
