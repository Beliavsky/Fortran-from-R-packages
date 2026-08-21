program test_rng
    use bgfd
    implicit none
    integer, parameter :: n = 20000
    real(dp) :: x(n), u(n), mean_u
    integer :: failures, i

    failures = 0
    call r_bell_ew(n, 0.9_dp, 1.6_dp, 1.3_dp, 0.7_dp, x)
    if (any(x < 0.0_dp)) failures = failures + 1
    do i = 1, n
        u(i) = p_bell_ew(x(i),0.9_dp,1.6_dp,1.3_dp,0.7_dp)
    end do
    mean_u = sum(u)/real(n,dp)
    if (abs(mean_u-0.5_dp) > 0.015_dp) then
        print '(a,f12.6)', 'Bell-EW PIT mean FAIL: ', mean_u
        failures = failures + 1
    end if

    call r_cbell_l(n, 1.5_dp, 2.3_dp, 0.7_dp, x)
    if (any(x < 0.0_dp)) failures = failures + 1
    do i = 1, n
        u(i) = p_cbell_l(x(i),1.5_dp,2.3_dp,0.7_dp)
    end do
    mean_u = sum(u)/real(n,dp)
    if (abs(mean_u-0.5_dp) > 0.015_dp) then
        print '(a,f12.6)', 'CBell-L PIT mean FAIL: ', mean_u
        failures = failures + 1
    end if

    if (failures /= 0) error stop 1
    print '(a)', 'test_rng: PASS'
end program test_rng
