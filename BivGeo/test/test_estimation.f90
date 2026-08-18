program test_estimation
    use bivgeo, only : dp, bivgeo_params, make_bivgeo_params, bivgeo_seed, rbivgeo2, mombivgeo
    implicit none

    integer, parameter :: n = 200000
    type(bivgeo_params) :: theta, estimate
    integer, allocatable :: z(:, :)
    logical :: ok

    allocate(z(n, 2))
    theta = make_bivgeo_params(0.45_dp, 0.60_dp, 0.72_dp)
    call bivgeo_seed(24680)
    call rbivgeo2(n, theta, z, ok)
    call check(ok, 1)
    call mombivgeo(z(:, 1), z(:, 2), estimate, ok)
    call check(ok, 2)
    call check(abs(estimate%theta1 - theta%theta1) < 0.025_dp, 3)
    call check(abs(estimate%theta2 - theta%theta2) < 0.025_dp, 4)
    call check(abs(estimate%theta3 - theta%theta3) < 0.025_dp, 5)

    print '(a)', 'test_estimation: PASS'

contains

    subroutine check(condition, code)
        logical, intent(in) :: condition
        integer, intent(in) :: code
        if (.not. condition) then
            print '(a,i0)', 'test_estimation: FAIL ', code
            error stop 1
        end if
    end subroutine check

end program test_estimation
