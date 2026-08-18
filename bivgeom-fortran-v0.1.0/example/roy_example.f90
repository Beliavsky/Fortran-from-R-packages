! Translation of bivgeom 1.0 computational code.
! Upstream DESCRIPTION declares License: GPL. See LICENSE.md.
program roy_example
    use, intrinsic :: iso_fortran_env, only : int64
    use bivgeom
    implicit none

    integer, allocatable :: sample(:, :)
    real(dp) :: theta(3)
    type(bivgeom_fit) :: fit

    theta = [0.5_dp, 0.7_dp, 0.9_dp]

    print '(a,f12.8)', 'P(X=2,Y=3) = ', dbivgeom_roy(2, 3, theta(1), theta(2), theta(3))
    print '(a,f12.8)', 'Corr(X,Y)   = ', corbivgeom_roy(theta(1), theta(2), theta(3))
    print '(a,f12.8)', 'P(X<=Y)     = ', relbivgeom_roy(theta(1), theta(2), theta(3))

    call seed_rng(int(12345, int64))
    call rbivgeom_roy(2000, theta(1), theta(2), theta(3), sample)
    fit = fit_bivgeom_ml(sample(:, 1), sample(:, 2))
    print '(a,3f10.5)', 'MLE theta = ', fit%theta

end program roy_example
