! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from R package good 1.0.2.

program example_good
    use good, only : dp, dgood, pgood, qgood, goodmean, good_glm_fit, glm_good
    implicit none

    integer :: y(10) = [0, 0, 0, 1, 1, 1, 1, 2, 2, 3]
    real(dp) :: x(10, 1)
    type(good_glm_fit) :: fit

    print '(a,f12.8)', 'P(X=4), z=0.6, s=-3: ', dgood(4, 0.6_dp, -3.0_dp)
    print '(a,f12.8)', 'P(X<=4):                 ', pgood(4.0_dp, 0.6_dp, -3.0_dp)
    print '(a,f12.3)',  'median:                  ', qgood(0.5_dp, 0.6_dp, -3.0_dp)
    print '(a,f12.8)', 'mean:                    ', goodmean(0.6_dp, -3.0_dp)

    x = 1.0_dp
    call glm_good(y, x, 'log', fit)
    print '(a,l1)', 'GLM converged: ', fit%converged
    print '(a,*(f12.6,1x))', 'coefficients [s, beta0]: ', fit%coefficients
    print '(a,f12.6)', 'log likelihood: ', fit%loglik
end program example_good
