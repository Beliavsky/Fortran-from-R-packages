! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

program test_simulation
    use dirmult_fortran
    implicit none
    real(dp), allocatable :: d(:,:)
    real(dp) :: alpha(3), p(3), means(3)
    integer :: info, fails, i
    type(sim_pop_result_type) :: sp
    type(null_test_result_type) :: nt
    type(dirmult_fit_type) :: fit
    type(mom_result_type) :: mr
    integer :: x(9,3)

    fails = 0
    alpha = [2.0_dp, 3.0_dp, 5.0_dp]
    call seed_rng(12345)
    call random_dirichlet(4000, alpha, d, info)
    call check_true('dirichlet info', info == 0, fails)
    do i = 1, size(d,1)
        if (abs(sum(d(i,:))-1.0_dp) > 2.0e-14_dp) fails = fails + 1
    end do
    means = sum(d, dim=1) / real(size(d,1),dp)
    call check_close('dir mean 1', means(1), 0.2_dp, 0.015_dp, fails)
    call check_close('dir mean 2', means(2), 0.3_dp, 0.015_dp, fails)
    call check_close('dir mean 3', means(3), 0.5_dp, 0.015_dp, fails)

    p = [0.2_dp, 0.3_dp, 0.5_dp]
    call sim_pop_equal_n(40, 3, 100, 0.1_dp, sp, pi=p, seed=777)
    call check_true('simPop info', sp%info == 0, fails)
    call check_true('simPop shape', size(sp%data,1) == 40 .and. size(sp%data,2) == 3, fails)
    call check_true('simPop totals', all(sum(sp%data,dim=2) == 100), fails)
    call check_true('simPop nonnegative', all(sp%data >= 0), fails)
    call check_close('simPop pi 2', sp%pi(2), 0.3_dp, 1.0e-15_dp, fails)

    x = reshape([ &
        25,22,3,4,2,3,20,5,5, &
        3,6,24,23,4,5,5,20,5, &
        2,2,3,3,24,22,5,5,20], shape(x))
    call null_test(x, nt, m=3, prec=8, store_data=.false., seed=2026)
    call check_true('null length', size(nt%mle_theta) == 4, fails)
    call check_true('null no stored cubes', .not. allocated(nt%simulated), fails)
    call fit_dirmult(x, fit, epsilon=1.0e-8_dp, trace=.false.)
    call weir_mom(x, mr)
    call check_close('null observed theta', nt%mle_theta(4), fit%theta, 1.0e-12_dp, fails)
    call check_close('null observed dm ll', nt%dm_loglik(4), fit%loglik, 1.0e-12_dp, fails)
    call check_close('null observed mom', nt%mom(4), mr%theta, 1.0e-12_dp, fails)
    call check_close('null observed mn', nt%mn_loglik(4), multinomial_loglik(x), 1.0e-12_dp, fails)

    if (fails /= 0) then
        write(*,'(a,i0)') 'test_simulation: FAIL ', fails
        error stop 1
    end if
    print *, 'test_simulation: PASS'

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

end program test_simulation
