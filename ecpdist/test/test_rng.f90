program test_rng
    use ecpdist, only: dp, recp, pecp, ecp_integral_result, ecp_kmoment
    implicit none
    integer, parameter :: n = 50000
    real(dp) :: x(n), ubar, xbar
    type(ecp_integral_result) :: m1
    integer :: seed_size, i, status, fails
    integer, allocatable :: seed(:)

    fails = 0
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    do i = 1, seed_size
        seed(i) = 104729 + 37*i
    end do
    call random_seed(put=seed)

    call recp(n, 1.0_dp, 1.0_dp, 1.0_dp, x, status)
    if (status /= 0) then
        print '(a,i0)', 'rng status: ', status
        fails = fails + 1
    end if
    if (minval(x) < 0.0_dp) then
        print '(a)', 'rng generated negative value'
        fails = fails + 1
    end if

    m1 = ecp_kmoment(1, 1.0_dp, 1.0_dp, 1.0_dp)
    xbar = sum(x)/real(n, dp)
    ubar = sum(pecp(x, 1.0_dp, 1.0_dp, 1.0_dp))/real(n, dp)
    if (abs(xbar - m1%estimate) > 0.015_dp) then
        print '(a,2f14.7)', 'rng mean mismatch: ', xbar, m1%estimate
        fails = fails + 1
    end if
    if (abs(ubar - 0.5_dp) > 0.006_dp) then
        print '(a,f14.7)', 'PIT mean mismatch: ', ubar
        fails = fails + 1
    end if

    if (fails /= 0) error stop 1
    print '(a)', 'test_rng: PASS'
end program test_rng
