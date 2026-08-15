program basic
    use delaporte, only : dp, ddelap, pdelap, qdelap, rdelap, momdelap, &
        seed_delaporte
    implicit none

    real(dp) :: x(10000), params(3)
    integer :: status

    print '(a,f0.12)', 'P(X=4) = ', ddelap(4.0_dp, 1.0_dp, 4.0_dp, 2.0_dp)
    print '(a,f0.12)', 'P(X<=4) = ', pdelap(4.0_dp, 1.0_dp, 4.0_dp, 2.0_dp)
    print '(a,f0.0)', '40th percentile = ', qdelap(0.4_dp, 1.0_dp, 4.0_dp, 2.0_dp)

    call seed_delaporte(4175)
    call rdelap(size(x), 10.0_dp, 2.0_dp, 10.0_dp, x)
    print '(a,f0.6)', 'sample mean = ', sum(x) / real(size(x), dp)

    call momdelap([5.0_dp, 7.0_dp, 9.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, &
        11.0_dp, 13.0_dp, 17.0_dp, 24.0_dp], params, 2, status)
    print '(a,3(1x,f0.8))', 'MoM alpha beta lambda =', params
end program basic
