program test_lmoments
    use lmoments, only: dp, lmoments_sample, lcoefs_sample, hosking_lmoments
    implicit none
    real(dp), parameter :: x(10) = [3.2_dp, -1.1_dp, 0.4_dp, 5.7_dp, 2.2_dp, &
        1.8_dp, 9.1_dp, -3.0_dp, 4.4_dp, 0.0_dp]
    real(dp), parameter :: ref(6) = [2.27_dp, 2.081111111111111_dp, &
        0.243333333333333_dp, 0.398571428571429_dp, 0.092222222222222_dp, 0.25_dp]
    real(dp) :: lmom(6), h(6), lc(6)
    integer :: info

    call lmoments_sample(x, lmom, info)
    call assert_true(info == 0, 'lmoments info')
    call assert_vec(lmom, ref, 2.0e-13_dp, 'lmoments reference')

    call hosking_lmoments(x, h, info)
    call assert_true(info == 0, 'hosking info')
    call assert_vec(h, ref, 2.0e-13_dp, 'Hosking reference path')

    call lcoefs_sample(x, lc, info)
    call assert_true(info == 0, 'lcoefs info')
    call assert_close(lc(1), ref(1), 2.0e-13_dp, 'L1')
    call assert_close(lc(2), ref(2), 2.0e-13_dp, 'L2')
    call assert_close(lc(3), ref(3) / ref(2), 2.0e-13_dp, 'tau3')
    call assert_close(lc(6), ref(6) / ref(2), 2.0e-13_dp, 'tau6')

    print '(a)', 'test_lmoments: PASS'
contains
    subroutine assert_vec(a, b, tol, label)
        real(dp), intent(in) :: a(:), b(:), tol
        character(*), intent(in) :: label
        if (maxval(abs(a - b)) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
    subroutine assert_close(a, b, tol, label)
        real(dp), intent(in) :: a, b, tol
        character(*), intent(in) :: label
        if (abs(a - b) > tol) then
            print *, label, a, b
            error stop 1
        end if
    end subroutine
    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, label
            error stop 1
        end if
    end subroutine
end program test_lmoments
