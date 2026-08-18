! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
program test_rng
    use, intrinsic :: iso_fortran_env, only : int64
    use bivgeom
    implicit none

    integer, allocatable :: z(:, :)
    integer :: i, n
    real(dp) :: mx, my, sx, sy, sxy, corr

    call seed_rng(int(12345, int64))
    n = 50000
    call rbivgeom_roy(n, 0.5_dp, 0.7_dp, 0.9_dp, z)
    if (size(z, 1) /= n .or. size(z, 2) /= 2) error stop 1
    if (any(z < 0)) error stop 1

    mx = sum(real(z(:, 1), dp)) / real(n, dp)
    my = sum(real(z(:, 2), dp)) / real(n, dp)
    if (abs(mx - 1.0_dp) > 0.035_dp) error stop 1
    if (abs(my - 7.0_dp / 3.0_dp) > 0.055_dp) error stop 1

    sx = 0.0_dp
    sy = 0.0_dp
    sxy = 0.0_dp
    do i = 1, n
        sx = sx + (real(z(i, 1), dp) - mx)**2
        sy = sy + (real(z(i, 2), dp) - my)**2
        sxy = sxy + (real(z(i, 1), dp) - mx) * (real(z(i, 2), dp) - my)
    end do
    corr = sxy / sqrt(sx * sy)
    if (abs(corr - corbivgeom_roy(0.5_dp, 0.7_dp, 0.9_dp)) > 0.025_dp) error stop 1

    print *, 'test_rng: PASS'

end program test_rng
