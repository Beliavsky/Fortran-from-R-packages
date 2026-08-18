! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
program test_estimation
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use bivgeom
    implicit none

    integer, allocatable :: z(:, :)
    real(dp) :: theta(3), truth(3), nll_truth
    type(bivgeom_fit) :: fit

    truth = [0.5_dp, 0.7_dp, 0.9_dp]
    call seed_rng(int(12345, int64))
    call rbivgeom_roy(10000, truth(1), truth(2), truth(3), z)

    theta = estbivgeom_roy(z(:, 1), z(:, 2), 'MMP')
    if (maxval(abs(theta(1:2) - truth(1:2))) > 0.025_dp) error stop 1
    if (abs(theta(3) - truth(3)) > 0.05_dp) error stop 1

    theta = estbivgeom_roy(z(:, 1), z(:, 2), 'LS')
    if (maxval(abs(theta(1:2) - truth(1:2))) > 0.035_dp) error stop 1
    if (abs(theta(3) - truth(3)) > 0.06_dp) error stop 1

    theta = estbivgeom_roy(z(:, 1), z(:, 2), 'MM4')
    if (.not. all(ieee_is_finite(theta))) error stop 1

    fit = fit_bivgeom_ml(z(1:2000, 1), z(1:2000, 2))
    if (.not. fit%converged) error stop 1
    if (.not. feasible_roy(fit%theta(1), fit%theta(2), fit%theta(3))) error stop 1
    if (maxval(abs(fit%theta - truth)) > 0.08_dp) error stop 1
    nll_truth = negative_loglik_roy(truth, z(1:2000, 1), z(1:2000, 2))
    if (fit%nll > nll_truth + 1.0e-8_dp) error stop 1

    print *, 'test_estimation: PASS'

end program test_estimation
