! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program tensor_example
    use mcmcglmm, only : dp, symmetrizer_matrix, normal_moment_matrix
    implicit none

    integer :: info
    real(dp), allocatable :: moments(:, :)
    real(dp), allocatable :: symmetrizer(:, :)

    call symmetrizer_matrix(2, 2, symmetrizer, info)
    if (info /= 0) error stop 'KPPM symmetrizer failed'
    call normal_moment_matrix(reshape([2.0_dp], [1, 1]), 4, moments, info)
    if (info /= 0) error stop 'normal moment tensor failed'

    print '(a,f8.4)', 'trace of KPPM(2,2): ', symmetrizer(1, 1) + symmetrizer(2, 2) + &
        symmetrizer(3, 3) + symmetrizer(4, 4)
    print '(a,f8.4)', 'E[X^4] for N(0,2): ', moments(1, 1)
end program tensor_example
