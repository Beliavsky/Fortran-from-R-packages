! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
program test_dependence
    use bivgeom
    implicit none

    real(dp), parameter :: t1 = 0.5_dp, t2 = 0.7_dp, t3 = 0.9_dp
    real(dp) :: cdf, mean_y, corr, rel, rel_direct, exy, mx, my, vx, vy, corr_direct
    real(dp) :: prob, prev
    integer :: x, y

    cdf = fyxbivgeom_roy(3.0_dp, t1, t2, t3, 2)
    mean_y = eyxbivgeom_roy(t1, t2, t3, 2)
    call check_close(cdf, 0.8611009774670882_dp, 2.0e-15_dp, 'conditional cdf')
    call check_close(mean_y, 1.576871072971574_dp, 3.0e-15_dp, 'conditional mean')

    mean_y = 0.0_dp
    prev = 0.0_dp
    do y = 0, 100
        cdf = fyxbivgeom_roy(real(y, dp), t1, t2, t3, 2)
        prob = cdf - prev
        mean_y = mean_y + real(y, dp) * prob
        prev = cdf
    end do
    call check_close(mean_y, eyxbivgeom_roy(t1, t2, t3, 2), 1.0e-12_dp, &
        'conditional mean from cdf')

    corr = corbivgeom_roy(t1, t2, t3)
    rel = relbivgeom_roy(t1, t2, t3)
    call check_close(corr, -0.23658256091460922_dp, 2.0e-15_dp, 'correlation reference')
    call check_close(rel, 0.7337952065042417_dp, 2.0e-15_dp, 'reliability reference')

    rel_direct = 0.0_dp
    exy = 0.0_dp
    do x = 0, 80
        do y = 0, 80
            prob = dbivgeom_roy(x, y, t1, t2, t3)
            if (x <= y) rel_direct = rel_direct + prob
            exy = exy + real(x * y, dp) * prob
        end do
    end do
    call check_close(rel_direct, rel, 5.0e-13_dp, 'reliability direct sum')

    mx = t1 / (1.0_dp - t1)
    my = t2 / (1.0_dp - t2)
    vx = t1 / (1.0_dp - t1)**2
    vy = t2 / (1.0_dp - t2)**2
    corr_direct = (exy - mx * my) / sqrt(vx * vy)
    call check_close(corr_direct, corr, 1.0e-12_dp, 'correlation direct sum')

    call check_close(lambda1_roy(2, t1, t3), 0.595_dp, 1.0e-15_dp, 'lambda1')
    call check_close(lambda2_roy(1, t2, t3), 0.37_dp, 1.0e-15_dp, 'lambda2')

    print *, 'test_dependence: PASS'

contains

    subroutine check_close(actual, expected, tol, label)
        real(dp), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print *, trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_close

end program test_dependence
