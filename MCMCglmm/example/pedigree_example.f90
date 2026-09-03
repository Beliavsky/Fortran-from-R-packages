! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program pedigree_example
    use mcmcglmm, only : dp, rng_state, rng_seed, pedigree_inverse, breeding_values_pedigree
    implicit none

    integer :: dam(4)
    integer :: info
    integer :: sire(4)
    real(dp), allocatable :: a_inverse(:, :)
    real(dp), allocatable :: breeding_values(:, :)
    real(dp), allocatable :: inbreeding(:)
    real(dp), allocatable :: mendelian_variance(:)
    real(dp) :: genetic_covariance(2, 2)
    type(rng_state) :: state

    dam = [0, 0, 1, 3]
    sire = [0, 0, 2, 2]
    call pedigree_inverse(dam, sire, a_inverse, inbreeding, mendelian_variance, info)
    if (info /= 0) error stop 'pedigree inverse failed'

    genetic_covariance = reshape([1.0_dp, 0.3_dp, 0.3_dp, 2.0_dp], [2, 2])
    call rng_seed(state, 20260831_8)
    call breeding_values_pedigree(state, dam, sire, genetic_covariance, breeding_values, info)
    if (info /= 0) error stop 'breeding-value simulation failed'

    print '(a,f8.4)', 'inbreeding of individual 4: ', inbreeding(4)
    print '(a,2f10.5)', 'breeding values of individual 4: ', breeding_values(4, :)
end program pedigree_example
