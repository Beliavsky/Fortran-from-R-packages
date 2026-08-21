! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

program basic_fit
    use dirmult_fortran
    implicit none
    integer :: x(9,3)
    integer :: i
    type(dirmult_fit_type) :: fit
    type(dirmult_summary_type) :: summary

    x = reshape([ &
        25,22,3,4,2,3,20,5,5, &
        3,6,24,23,4,5,5,20,5, &
        2,2,3,3,24,22,5,5,20], shape(x))

    call fit_dirmult(x, fit, epsilon=1.0e-10_dp, trace=.false.)
    call summarize_dirmult(x, fit, summary)

    write(*,'(a,f12.8)') 'theta = ', fit%theta
    write(*,'(a,es18.9)') 'loglik = ', fit%loglik
    do i = 1, size(fit%pi)
        write(*,'(a,i0,a,f12.8,a,f12.8)') 'pi(', i, ') = ', fit%pi(i), &
            '  se = ', summary%se_mle(i)
    end do
end program basic_fit
