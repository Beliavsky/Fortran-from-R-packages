program test_rng_gof
    use bivpois, only : dp, seed_rng, rbp, bp_gof_result, bp_gof, bp_gof2
    implicit none
    integer, parameter :: n = 40000
    real(dp), parameter :: lambda(3) = [3.0_dp, 5.0_dp, 2.0_dp]
    integer, allocatable :: x(:, :)
    real(dp) :: m1, m2, cov
    real(dp), parameter :: lambda_large(3) = [40.0_dp, 50.0_dp, 10.0_dp]
    type(bp_gof_result) :: g1, g2
    integer, parameter :: a(20) = [5,6,4,7,3,5,8,4,6,5,7,3,5,4,6,8,5,4,7,6]
    integer, parameter :: b(20) = [7,8,6,9,5,8,9,6,8,7,9,5,7,6,8,10,7,6,9,8]

    allocate(x(n, 2))
    call seed_rng(13579)
    call rbp(n, lambda, x)
    m1 = sum(real(x(:, 1), dp)) / real(n, dp)
    m2 = sum(real(x(:, 2), dp)) / real(n, dp)
    cov = sum((real(x(:, 1), dp) - m1) * (real(x(:, 2), dp) - m2)) / real(n - 1, dp)
    if (abs(m1 - 5.0_dp) > 0.06_dp) error stop "rbp first marginal mean"
    if (abs(m2 - 7.0_dp) > 0.07_dp) error stop "rbp second marginal mean"
    if (abs(cov - 2.0_dp) > 0.12_dp) error stop "rbp covariance"


    call seed_rng(97531)
    call rbp(n, lambda_large, x)
    m1 = sum(real(x(:, 1), dp)) / real(n, dp)
    m2 = sum(real(x(:, 2), dp)) / real(n, dp)
    cov = sum((real(x(:, 1), dp) - m1) * (real(x(:, 2), dp) - m2)) / real(n - 1, dp)
    if (abs(m1 - 50.0_dp) > 0.15_dp) error stop "large-rbp first marginal mean"
    if (abs(m2 - 60.0_dp) > 0.17_dp) error stop "large-rbp second marginal mean"
    if (abs(cov - 10.0_dp) > 0.45_dp) error stop "large-rbp covariance"

    call seed_rng(24680)
    g1 = bp_gof(a, b, 39)
    if (g1%pvalue <= 0.0_dp .or. g1%pvalue > 1.0_dp) error stop "bp_gof pvalue"
    if (g1%replicates /= 39) error stop "bp_gof replicates"

    call seed_rng(24680)
    g2 = bp_gof2(a, b, 39)
    if (abs(g1%pvalue - g2%pvalue) > 0.0_dp) error stop "bp_gof2 parity"
    if (maxval(abs(g1%lambda - g2%lambda)) > 0.0_dp) error stop "bp_gof2 fit parity"

    print '(a)', 'test_rng_gof: PASS'
end program test_rng_gof
