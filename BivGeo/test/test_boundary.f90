program test_boundary
    use bivgeo, only : dp, bivgeo_params, make_bivgeo_params, dbivgeo2, covbivgeo, corbivgeo, bivgeo_seed, rbivgeo2
    implicit none

    integer, parameter :: n = 100000
    type(bivgeo_params) :: theta
    integer, allocatable :: z(:, :)
    real(dp) :: expected, mx, my, cxy, nr
    logical :: ok

    theta = make_bivgeo_params(0.3_dp, 0.6_dp, 1.0_dp)
    expected = (1.0_dp - 0.3_dp) * 0.3_dp**1 * (1.0_dp - 0.6_dp) * 0.6_dp**2
    call check(abs(dbivgeo2(2, 3, theta) - expected) < 1.0e-14_dp, 1)
    call check(abs(covbivgeo(theta)) < 1.0e-14_dp, 2)
    call check(abs(corbivgeo(theta)) < 1.0e-14_dp, 3)

    allocate(z(n, 2))
    call bivgeo_seed(777)
    call rbivgeo2(n, theta, z, ok)
    call check(ok, 4)
    nr = real(n, dp)
    mx = sum(real(z(:, 1), dp)) / nr
    my = sum(real(z(:, 2), dp)) / nr
    cxy = sum((real(z(:, 1), dp) - mx) * (real(z(:, 2), dp) - my)) / (nr - 1.0_dp)
    call check(abs(cxy) < 0.025_dp, 5)

    print '(a)', 'test_boundary: PASS'

contains

    subroutine check(condition, code)
        logical, intent(in) :: condition
        integer, intent(in) :: code
        if (.not. condition) then
            print '(a,i0)', 'test_boundary: FAIL ', code
            error stop 1
        end if
    end subroutine check

end program test_boundary
