! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
program test_distribution
    use bivgeom
    implicit none

    real(dp), parameter :: t1 = 0.5_dp, t2 = 0.7_dp, t3 = 0.9_dp
    real(dp) :: p, f, s, total, cell
    integer :: x, y

    call check_close(dbivgeom_roy(2, 3, t1, t2, t3), &
        0.011598034467511013_dp, 2.0e-15_dp, 'pmf reference')
    call check_close(fbivgeom_roy(2.0_dp, 3.0_dp, t1, t2, t3), &
        0.6433764164636361_dp, 2.0e-15_dp, 'cdf reference')
    call check_close(sbivgeom_roy(2.0_dp, 3.0_dp, t1, t2, t3), &
        0.04557106575_dp, 2.0e-15_dp, 'survival reference')

    p = dbivgeom_roy(2, 3, t1, t2, t3)
    f = fbivgeom_roy(2.0_dp, 3.0_dp, t1, t2, t3)
    s = sbivgeom_roy(2.0_dp, 3.0_dp, t1, t2, t3)
    if (p <= 0.0_dp .or. f <= 0.0_dp .or. s <= 0.0_dp) error stop 1

    cell = fbivgeom_roy(2.0_dp, 3.0_dp, t1, t2, t3) - &
        fbivgeom_roy(1.0_dp, 3.0_dp, t1, t2, t3) - &
        fbivgeom_roy(2.0_dp, 2.0_dp, t1, t2, t3) + &
        fbivgeom_roy(1.0_dp, 2.0_dp, t1, t2, t3)
    call check_close(cell, p, 5.0e-15_dp, 'cdf finite difference')

    total = 0.0_dp
    do x = 0, 60
        do y = 0, 80
            total = total + dbivgeom_roy(x, y, t1, t2, t3)
        end do
    end do
    call check_close(total, 1.0_dp, 2.0e-9_dp, 'pmf normalization')

    if (.not. feasible_roy(t1, t2, t3)) error stop 1
    if (feasible_roy(0.5_dp, 1.2_dp, 0.9_dp)) error stop 1
    if (feasible_roy(0.9_dp, 0.9_dp, 0.5_dp)) error stop 1

    print *, 'test_distribution: PASS'

contains

    subroutine check_close(actual, expected, tol, label)
        real(dp), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print *, trim(label), actual, expected
            error stop 1
        end if
    end subroutine check_close

end program test_distribution
