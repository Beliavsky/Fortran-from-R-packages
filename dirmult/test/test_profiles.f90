! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

program test_profiles
    use dirmult_fortran
    implicit none
    integer :: x1(9,3), x2(8,3), fails
    type(profile_fit_type) :: pf
    type(profile_grid_type) :: grid
    type(count_table_type) :: tabs(2)
    type(equal_theta_fit_type) :: eq

    fails = 0
    x1 = reshape([ &
        25,22,3,4,2,3,20,5,5, &
        3,6,24,23,4,5,5,20,5, &
        2,2,3,3,24,22,5,5,20], shape(x1))
    x2 = reshape([ &
        18,20,4,5,3,4,15,6, &
        8,7,20,19,5,6,10,16, &
        4,3,6,6,22,20,5,8], shape(x2))

    call estimate_profile_loglik(x1, 0.2_dp, pf, epsilon=1.0e-10_dp, trace=.false.)
    call check_true('profile converged', pf%converged, fails)
    call check_close('profile ll', pf%loglik, -227.33355744330612_dp, 5.0e-9_dp, fails)
    call check_close('profile pi1', pf%pi(1), 0.3254196899158144_dp, 2.0e-10_dp, fails)
    call check_close('profile theta', pf%theta, 0.2_dp, 2.0e-10_dp, fails)

    call grid_profile(x1, 0.2_dp, -0.05_dp, 0.05_dp, 3, grid, epsilon=1.0e-9_dp)
    call check_true('grid length', size(grid%theta) == 3, fails)
    call check_close('grid center theta', grid%theta(2), 0.2_dp, 1.0e-15_dp, fails)
    call check_close('grid center ll', grid%loglik(2), -227.33355744330612_dp, 2.0e-8_dp, fails)

    allocate(tabs(1)%x(9,3), tabs(2)%x(8,3))
    tabs(1)%x = x1
    tabs(2)%x = x2
    call fit_equal_theta(tabs, 0.25_dp, eq, epsilon=1.0e-10_dp, trace=.false.)
    call check_true('equal theta converged', eq%converged, fails)
    call check_close('equal loglik', eq%loglik, -464.58656986221547_dp, 2.0e-7_dp, fails)
    call check_close('equal theta 1', eq%table(1)%theta, 0.1942962982040319_dp, 2.0e-10_dp, fails)
    call check_close('equal theta 2', eq%table(2)%theta, 0.1942962982040319_dp, 2.0e-10_dp, fails)
    call check_close('equal pi 1.1', eq%table(1)%pi(1), 0.325350432_dp, 5.0e-8_dp, fails)
    call check_close('equal pi 2.2', eq%table(2)%pi(2), 0.38417484_dp, 3.0e-7_dp, fails)

    if (fails /= 0) then
        write(*,'(a,i0)') 'test_profiles: FAIL ', fails
        error stop 1
    end if
    print *, 'test_profiles: PASS'

contains

    subroutine check_close(name, got, want, tol, nf)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        integer, intent(inout) :: nf
        if (abs(got-want) > tol) then
            write(*,'(a,2(1x,es24.15))') trim(name)//' got/want:', got, want
            nf = nf + 1
        end if
    end subroutine check_close

    subroutine check_true(name, ok, nf)
        character(len=*), intent(in) :: name
        logical, intent(in) :: ok
        integer, intent(inout) :: nf
        if (.not. ok) then
            write(*,'(a)') trim(name)//': failed'
            nf = nf + 1
        end if
    end subroutine check_true

end program test_profiles
