! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

program test_distribution
    use good, only : dp, dgood, pgood, qgood, goodmean, good_moments, rgood
    implicit none

    real(dp), parameter :: tol = 2.0e-10_dp
    real(dp) :: z, s, p, mu, var, exact
    integer :: status, sample(20000)
    real(dp) :: sample_mean

    z = 0.6_dp

    ! s = 0 is geometric on x=0,1,... with P(X=x)=(1-z) z^x.
    s = 0.0_dp
    exact = (1.0_dp - z) * z**4
    call assert_close(dgood(4, z, s), exact, tol, 'dgood s=0')
    call assert_close(pgood(3.0_dp, z, s), 1.0_dp - z**4, tol, 'pgood s=0')
    call assert_close(goodmean(z, s), z / (1.0_dp - z), tol, 'goodmean s=0')
    call good_moments(z, s, mu, var, status)
    call assert_close(mu, z / (1.0_dp - z), tol, 'moment mean s=0')
    call assert_close(var, z / (1.0_dp - z)**2, 2.0e-9_dp, 'moment variance s=0')

    ! s = -1 gives a negative-binomial form: (x+1)(1-z)^2 z^x.
    s = -1.0_dp
    exact = 5.0_dp * (1.0_dp - z)**2 * z**4
    call assert_close(dgood(4, z, s), exact, 2.0e-10_dp, 'dgood s=-1')
    call assert_close(goodmean(z, s), 2.0_dp * z / (1.0_dp - z), 2.0e-9_dp, 'goodmean s=-1')

    ! Quantile inversion consistency.
    s = -3.0_dp
    p = 0.73_dp
    exact = qgood(p, z, s, status=status)
    if (status < 0) error stop 'qgood returned error'
    if (pgood(exact, z, s) < p) error stop 'qgood lower CDF consistency failed'
    if (exact > 0.0_dp) then
        if (pgood(exact - 1.0_dp, z, s) >= p) error stop 'qgood upper CDF consistency failed'
    end if

    ! Exercise the extreme-s log-sum-exp path used instead of overflowing polylog.
    p = dgood(330, 0.6_dp, -170.0_dp, status)
    if (status < 0 .or. p <= 0.0_dp .or. p >= 1.0_dp) error stop 'extreme-s dgood failed'

    call rgood(size(sample), 0.4_dp, 0.0_dp, sample, status=status)
    if (status < 0) error stop 'rgood failed'
    sample_mean = sum(real(sample, dp)) / real(size(sample), dp)
    if (abs(sample_mean - 0.4_dp / 0.6_dp) > 0.06_dp) error stop 'rgood sample mean failed'

    print '(a)', 'test_distribution: PASS'

contains

    subroutine assert_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual, expected, atol
        character(len=*), intent(in) :: label
        if (abs(actual - expected) > atol) then
            print '(a,2es24.15)', trim(label)//' actual/expected: ', actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_distribution
