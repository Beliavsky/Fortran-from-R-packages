program test_distributions
    use rangen
    implicit none
    integer, parameter :: n = 80000
    real(dp), allocatable :: x(:)
    real(dp) :: m, v, target

    call seed_all(12345_i8)

    x = runif(n, 2.0_dp, 4.0_dp)
    call check_mean(x, 3.0_dp, 0.015_dp, "uniform mean")

    x = rbeta(n, 2.0_dp, 5.0_dp)
    call check_mean(x, 2.0_dp / 7.0_dp, 0.006_dp, "beta mean")

    x = rexp(n, 2.0_dp)
    call check_mean(x, 0.5_dp, 0.01_dp, "exponential mean")

    x = rchisq(n, 4.0_dp)
    call check_mean(x, 4.0_dp, 0.04_dp, "chisq mean")

    x = rgamma(n, 3.0_dp, 2.0_dp)
    call check_mean(x, 1.5_dp, 0.02_dp, "gamma mean")

    x = rgamma(n, 0.4_dp, 1.7_dp)
    call check_mean(x, 0.4_dp / 1.7_dp, 0.008_dp, "gamma shape below one mean")

    x = rbeta(n, 0.5_dp, 0.8_dp)
    call check_mean(x, 0.5_dp / 1.3_dp, 0.008_dp, "beta subunit shapes mean")

    x = rgeom(n, 0.25_dp)
    call check_mean(x, 3.0_dp, 0.05_dp, "geometric mean")

    x = rnorm(n, 1.0_dp, 2.0_dp)
    m = sum(x) / real(n, dp)
    v = sum((x - m) ** 2) / real(n - 1, dp)
    if (abs(m - 1.0_dp) > 0.025_dp) error stop "normal mean"
    if (abs(v - 4.0_dp) > 0.08_dp) error stop "normal variance"

    x = rt(n, 8.0_dp, 0.0_dp)
    call check_mean(x, 0.0_dp, 0.025_dp, "t mean")

    x = rpareto(n, 3.0_dp, 2.0_dp)
    call check_mean(x, 3.0_dp, 0.05_dp, "pareto mean")

    x = rfrechet(n, 4.0_dp, 1.0_dp, 2.0_dp)
    target = 1.0_dp + 2.0_dp * gamma(0.75_dp)
    call check_mean(x, target, 0.05_dp, "frechet mean")

    x = rlaplace(n, 1.0_dp, 2.0_dp)
    call check_mean(x, 1.0_dp, 0.035_dp, "laplace mean")

    x = rgumbel(n, 1.0_dp, 2.0_dp)
    call check_mean(x, 1.0_dp + 2.0_dp * euler_gamma, 0.035_dp, "gumbel mean")

    x = rarcsine(n, -1.0_dp, 3.0_dp)
    call check_mean(x, 1.0_dp, 0.02_dp, "arcsine mean")

    x = rcauchy(n, 2.0_dp, 1.0_dp)
    if (abs(real(count(x <= 2.0_dp), dp) / real(n, dp) - 0.5_dp) > 0.01_dp) then
        error stop "cauchy median probability"
    end if

    print *, "test_distributions: PASS"

contains

    subroutine check_mean(a, expected, tol, label)
        real(dp), intent(in) :: a(:), expected, tol
        character(*), intent(in) :: label
        real(dp) :: amean
        amean = sum(a) / real(size(a), dp)
        if (abs(amean - expected) > tol) then
            print *, trim(label), amean, expected
            error stop "mean check failed"
        end if
    end subroutine check_mean


end program test_distributions
