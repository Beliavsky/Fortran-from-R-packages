program test_skewt
    use skewt, only : dp, dskt, pskt, qskt, rskt
    use skewt_special, only : student_t_pdf, student_t_cdf
    implicit none

    integer, parameter :: n = 40000
    real(dp) :: probs(9), q(9), xgrid(9), vals(9)
    real(dp), allocatable :: r(:)
    real(dp) :: m, pneg, expected_pneg
    integer :: i, failures

    failures = 0

    ! gamma = 1 must reduce exactly to the ordinary Student-t law.
    xgrid = [-3.0_dp, -1.0_dp, -0.1_dp, 0.0_dp, 0.1_dp, 0.5_dp, &
        1.0_dp, 2.0_dp, 4.0_dp]
    vals = dskt(xgrid, 5.0_dp, 1.0_dp)
    do i = 1, size(xgrid)
        call check_close(vals(i), student_t_pdf(xgrid(i), 5.0_dp), &
            2.0e-13_dp, "symmetric density", failures)
        call check_close(pskt(xgrid(i), 5.0_dp, 1.0_dp), &
            student_t_cdf(xgrid(i), 5.0_dp), 2.0e-12_dp, &
            "symmetric cdf", failures)
    end do

    ! Values from the formulas/examples in the R package.
    call check_close(dskt(0.5_dp, 2.0_dp), 0.2962962962962963_dp, &
        5.0e-14_dp, "dskt example", failures)
    call check_close(dskt(0.01_dp, 2.0_dp, 2.0_dp), &
        0.282837409256623_dp, 2.0e-13_dp, "dskt skew example", failures)
    call check_close(pskt(1.25_dp, 2.0_dp, 2.0_dp), &
        0.5233808333817773_dp, 3.0e-12_dp, "pskt example", failures)
    call check_close(qskt(0.975_dp, 2.0_dp, 2.0_dp), &
        11.046797999046873_dp, 2.0e-9_dp, "qskt reference", failures)

    ! CDF/quantile inversion, including both sides of zero.
    probs = [1.0e-5_dp, 0.001_dp, 0.025_dp, 0.10_dp, 0.20_dp, &
        0.50_dp, 0.75_dp, 0.975_dp, 0.9999_dp]
    q = qskt(probs, 2.5_dp, 2.0_dp)
    do i = 1, size(probs)
        call check_close(pskt(q(i), 2.5_dp, 2.0_dp), probs(i), &
            2.0e-10_dp, "cdf quantile inversion", failures)
    end do

    ! P(X < 0) is the characteristic Fernandez-Steel mass split.
    expected_pneg = 1.0_dp / 5.0_dp
    call check_close(pskt(0.0_dp, 4.0_dp, 2.0_dp), expected_pneg, &
        2.0e-14_dp, "mass below zero", failures)

    allocate(r(n))
    call rskt(r, 7.0_dp, 2.0_dp)
    pneg = real(count(r < 0.0_dp), dp) / real(n, dp)
    call check_close(pneg, 0.2_dp, 0.012_dp, "rng negative mass", failures)

    ! Mean for df > 1: E|T_df| * (gamma^2 - gamma^-2)/(gamma + gamma^-1).
    ! E|T_df| = 2*sqrt(df)*Gamma((df+1)/2)/(sqrt(pi)*(df-1)*Gamma(df/2)).
    m = sum(r) / real(n, dp)
    call check_close(m, fs_mean(7.0_dp, 2.0_dp), 0.04_dp, &
        "rng mean", failures)

    if (failures == 0) then
        print '(a)', 'test_skewt: PASS'
    else
        print '(a,i0)', 'test_skewt: FAIL ', failures
        error stop 1
    end if

contains

    real(dp) function fs_mean(df, gamma) result(mu)
        real(dp), intent(in) :: df, gamma
        real(dp) :: eat
        eat = 2.0_dp * sqrt(df) * exp(log_gamma(0.5_dp * (df + 1.0_dp)) &
            - log_gamma(0.5_dp * df)) / (sqrt(acos(-1.0_dp)) * (df - 1.0_dp))
        mu = eat * (gamma * gamma - 1.0_dp / (gamma * gamma)) &
            / (gamma + 1.0_dp / gamma)
    end function fs_mean

    subroutine check_close(actual, expected, tol, label, failures)
        real(dp), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        if (abs(actual - expected) > tol) then
            print '(a,2(1x,es24.15),1x,a,es12.4)', trim(label), actual, expected, &
                'tol=', tol
            failures = failures + 1
        end if
    end subroutine check_close

end program test_skewt
