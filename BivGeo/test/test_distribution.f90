program test_distribution
    use bivgeo, only : dp, bivgeo_params, make_bivgeo_params, dbivgeo1, dbivgeo2, pbivgeo, sbivgeo
    implicit none

    type(bivgeo_params) :: theta
    real(dp) :: p1, p2, total, partial, tol
    integer :: x, y

    tol = 5.0e-13_dp
    theta = make_bivgeo_params(0.2_dp, 0.4_dp, 0.7_dp)

    p1 = dbivgeo1(1, 2, theta)
    p2 = dbivgeo2(1, 2, theta)
    call check(abs(p1 - 0.16128_dp) < tol, 1)
    call check(abs(p2 - 0.16128_dp) < tol, 2)
    call check(abs(dbivgeo1(1, 2, theta, .true.) - log(0.16128_dp)) < tol, 3)
    call check(abs(pbivgeo(1, 2, theta) - 0.79728_dp) < tol, 4)
    call check(abs(sbivgeo(1, 2, theta) - 0.01568_dp) < tol, 5)

    total = 0.0_dp
    do x = 1, 120
        do y = 1, 120
            p1 = dbivgeo1(x, y, theta)
            p2 = dbivgeo2(x, y, theta)
            call check(abs(p1 - p2) < 5.0e-14_dp, 6)
            total = total + p2
        end do
    end do
    call check(abs(total - 1.0_dp) < 2.0e-13_dp, 7)

    partial = 0.0_dp
    do x = 1, 3
        do y = 1, 4
            partial = partial + dbivgeo2(x, y, theta)
        end do
    end do
    call check(abs(partial - pbivgeo(3, 4, theta)) < 2.0e-13_dp, 8)

    print '(a)', 'test_distribution: PASS'

contains

    subroutine check(condition, code)
        logical, intent(in) :: condition
        integer, intent(in) :: code
        if (.not. condition) then
            print '(a,i0)', 'test_distribution: FAIL ', code
            error stop 1
        end if
    end subroutine check

end program test_distribution
