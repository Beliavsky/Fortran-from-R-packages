program test_rng
    use bivgeo, only : dp, bivgeo_params, make_bivgeo_params, mean_bivgeo, covbivgeo, bivgeo_seed, rbivgeo1, rbivgeo2
    implicit none

    integer, parameter :: n = 150000
    type(bivgeo_params) :: theta
    integer, allocatable :: z(:, :)
    real(dp) :: mu(2), mx, my, cxy
    logical :: ok

    allocate(z(n, 2))
    theta = make_bivgeo_params(0.5_dp, 0.5_dp, 0.7_dp)
    mu = mean_bivgeo(theta)

    call bivgeo_seed(12345)
    call rbivgeo2(n, theta, z, ok)
    call check(ok, 1)
    call sample_stats(z, mx, my, cxy)
    call check(abs(mx - mu(1)) < 0.015_dp, 2)
    call check(abs(my - mu(2)) < 0.015_dp, 3)
    call check(abs(cxy - covbivgeo(theta)) < 0.025_dp, 4)

    call bivgeo_seed(54321)
    call rbivgeo1(n, theta, z, ok)
    call check(ok, 5)
    call sample_stats(z, mx, my, cxy)
    call check(abs(mx - mu(1)) < 0.015_dp, 6)
    call check(abs(my - mu(2)) < 0.015_dp, 7)
    call check(abs(cxy - covbivgeo(theta)) < 0.025_dp, 8)

    print '(a)', 'test_rng: PASS'

contains

    subroutine sample_stats(a, mx, my, covxy)
        integer, intent(in) :: a(:, :)
        real(dp), intent(out) :: mx, my, covxy
        real(dp) :: nreal

        nreal = real(size(a, 1), dp)
        mx = sum(real(a(:, 1), dp)) / nreal
        my = sum(real(a(:, 2), dp)) / nreal
        covxy = sum((real(a(:, 1), dp) - mx) * (real(a(:, 2), dp) - my)) / (nreal - 1.0_dp)
    end subroutine sample_stats

    subroutine check(condition, code)
        logical, intent(in) :: condition
        integer, intent(in) :: code
        if (.not. condition) then
            print '(a,i0)', 'test_rng: FAIL ', code
            error stop 1
        end if
    end subroutine check

end program test_rng
