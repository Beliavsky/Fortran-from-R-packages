! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

program test_core
    use dirmult_fortran
    implicit none
    integer :: x(9,3), xz(10,4), fails
    real(dp) :: f(3,3), ef(3,3), s(3)
    type(mom_result_type) :: mr
    type(dirmult_fit_type) :: fit, fit_exp, fit_zero
    type(dirmult_summary_type) :: sm

    fails = 0
    x = reshape([ &
        25,22,3,4,2,3,20,5,5, &
        3,6,24,23,4,5,5,20,5, &
        2,2,3,3,24,22,5,5,20], shape(x))

    call weir_mom(x, mr)
    call check_close('weir theta', mr%theta, 0.39145920495712516_dp, 1.0e-12_dp, fails)
    call check_close('weir se', mr%se, 0.12528147955958720_dp, 1.0e-12_dp, fails)

    call fit_dirmult(x, fit, epsilon=1.0e-12_dp, trace=.false.)
    call check_true('fit converged', fit%converged, fails)
    call check_close('gamma 1', fit%gamma(1), 1.03582212_dp, 5.0e-8_dp, fails)
    call check_close('gamma 2', fit%gamma(2), 1.15481351_dp, 5.0e-8_dp, fails)
    call check_close('gamma 3', fit%gamma(3), 0.98728973_dp, 5.0e-8_dp, fails)
    call check_close('theta mle', fit%theta, 0.23935324699743565_dp, 2.0e-10_dp, fails)
    call check_close('loglik', fit%loglik, -227.06998093951688_dp, 2.0e-10_dp, fails)

    call fit_dirmult(x, fit_exp, epsilon=1.0e-10_dp, trace=.false., mode='exp')
    call check_true('expected fit converged', fit_exp%converged, fails)
    call check_close('expected theta', fit_exp%theta, fit%theta, 5.0e-8_dp, fails)

    xz = 0
    xz(1:9,1:3) = x
    call fit_dirmult(xz, fit_zero, epsilon=1.0e-12_dp, trace=.false.)
    call check_true('zero row/col removed', fit_zero%converged, fails)
    call check_close('zero filtered theta', fit_zero%theta, fit%theta, 2.0e-12_dp, fails)

    s = score_function(x, fit%gamma)
    call check_true('score near zero', maxval(abs(s)) < 2.0e-10_dp, fails)
    f = observed_fim(x, fit%gamma)
    call check_close('obs fim 11', f(1,1), -9.35352120_dp, 2.0e-7_dp, fails)
    call check_close('obs fim 12', f(1,2), 3.04809127_dp, 2.0e-7_dp, fails)
    ef = expected_fim(x, fit%gamma)
    call check_close('exp fim 11', ef(1,1), 8.77505946_dp, 3.0e-7_dp, fails)
    call check_close('exp fim 12', ef(1,2), -3.04809127_dp, 2.0e-7_dp, fails)

    call summarize_dirmult(x, fit, sm)
    call check_true('summary info', sm%info == 0, fails)
    call check_close('se mle pi1', sm%se_mle(1), 0.07350489_dp, 2.0e-8_dp, fails)
    call check_close('se mle pi2', sm%se_mle(2), 0.07606886_dp, 2.0e-8_dp, fails)
    call check_close('se mle pi3', sm%se_mle(3), 0.07226832_dp, 2.0e-8_dp, fails)
    call check_close('se mle theta', sm%se_mle(4), 0.05839951_dp, 2.0e-8_dp, fails)
    call check_close('mom pi1', sm%mom(1), 0.3296296296296296_dp, 1.0e-14_dp, fails)
    call check_close('mom se pi1', sm%se_mom(1), 0.31554382_dp, 2.0e-8_dp, fails)
    call check_close('mn loglik', multinomial_loglik(x), -296.39311984354765_dp, 2.0e-10_dp, fails)

    if (fails /= 0) then
        write(*,'(a,i0)') 'test_core: FAIL ', fails
        error stop 1
    end if
    print *, 'test_core: PASS'

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

end program test_core
