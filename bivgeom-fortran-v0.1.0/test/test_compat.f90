! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
program test_compat
    use, intrinsic :: iso_fortran_env, only : int64
    use bivgeom
    implicit none

    integer, allocatable :: z(:, :)
    real(dp) :: theta(3), a, b

    a = dbivgeomroy(2, 3, 0.5_dp, 0.7_dp, 0.9_dp)
    b = dbivgeom_roy(2, 3, 0.5_dp, 0.7_dp, 0.9_dp)
    if (abs(a - b) > 1.0e-15_dp) error stop 1
    if (abs(lambda1roy(7, 2, 0.5_dp, 0.7_dp, 0.9_dp) - 0.595_dp) > 1.0e-15_dp) error stop 1
    if (abs(lambda2roy(1, 9, 0.5_dp, 0.7_dp, 0.9_dp) - 0.37_dp) > 1.0e-15_dp) error stop 1

    call seed_rng(int(9, int64))
    call rbivgeomroy(500, 0.5_dp, 0.7_dp, 0.9_dp, z)
    theta = estbivgeomroy(z(:, 1), z(:, 2), 'MMP')
    if (any(theta <= 0.0_dp)) error stop 1
    if (s_n(1, 1, z(:, 1), z(:, 2)) < 0.0_dp) error stop 1
    if (minuslogroy(z(:, 1), z(:, 2), 0.5_dp, 0.7_dp, 0.9_dp) <= 0.0_dp) error stop 1

    print *, 'test_compat: PASS'

end program test_compat
