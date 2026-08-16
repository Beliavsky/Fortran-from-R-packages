program test_numerics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use sadists
    implicit none
    integer :: failures
    real(dp) :: k(6), m(6), c(6), p, q

    failures = 0

    call check(chisq_log_moment(10.0_dp, 3.0_dp, 1.0_dp), &
        2.5649493574615363_dp, 2.0e-13_dp, 'chisq moment 1', failures)
    call check(chisq_log_moment(10.0_dp, 3.0_dp, -1.0_dp), &
        -2.3518988421755163_dp, 2.0e-13_dp, 'chisq moment -1', failures)
    call check(digamma_dp(3.7_dp), 1.1671535393615113_dp, &
        2.0e-13_dp, 'digamma', failures)
    call check(polygamma_dp(1, 3.7_dp), 0.31003785767003833_dp, &
        2.0e-13_dp, 'trigamma', failures)
    call check(polygamma_dp(2, 3.7_dp), -0.09539530872855405_dp, &
        2.0e-13_dp, 'polygamma2', failures)

    m = [1.0_dp, 3.0_dp, 10.0_dp, 40.0_dp, 200.0_dp, 1200.0_dp]
    call moments_to_cumulants(m, c)
    call cumulants_to_moments(c, k)
    if (maxval(abs(k-m)) > 1.0e-11_dp) then
        print *, 'moment/cumulant round trip failed'
        failures = failures + 1
    end if

    k = [1.2_dp, 2.3_dp, 0.8_dp, 1.5_dp, 2.1_dp, 3.2_dp]
    call check(edgeworth_pdf(0.8_dp, k), 0.26808874089298024_dp, &
        2.0e-13_dp, 'edgeworth pdf', failures)
    call check(edgeworth_cdf(0.8_dp, k), 0.40654001802313927_dp, &
        2.0e-13_dp, 'edgeworth cdf', failures)
    call check(cornish_fisher_quantile(0.7_dp, k), &
        1.9332801310192802_dp, 2.0e-13_dp, 'cornish fisher', failures)
    call check(as269(0.2_dp, [0.5_dp, 0.2_dp, 0.1_dp, 0.05_dp]), &
        0.12218422438271608_dp, 2.0e-13_dp, 'AS269', failures)

    q = normal_quantile(-50.0_dp, lower_tail=.false., log_p=.true.)
    p = normal_cdf(q, lower_tail=.false.)
    if (.not. ieee_is_finite(q) .or. abs(log(p) + 50.0_dp) > 2.0e-6_dp) then
        print *, 'log upper-tail normal quantile failed', q, p
        failures = failures + 1
    end if

    if (failures == 0) then
        print *, 'test_numerics: PASS'
    else
        print *, 'test_numerics: FAIL', failures
        error stop 1
    end if

contains

    subroutine check(actual, expected, tol, name, failures)
        real(dp), intent(in) :: actual, expected, tol
        character(*), intent(in) :: name
        integer, intent(inout) :: failures
        if (abs(actual - expected) > tol * max(1.0_dp, abs(expected))) then
            print '(a,2es24.16)', trim(name)//': ', actual, expected
            failures = failures + 1
        end if
    end subroutine check

end program test_numerics
