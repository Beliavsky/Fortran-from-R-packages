! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program distributions_example
    use mcmcglmm, only : dp, rng_state, rng_seed, truncated_normal_sample, riw_mcmcglmm, pkk_probability
    implicit none

    integer :: info
    real(dp) :: draw
    real(dp), allocatable :: covariance(:, :)
    real(dp) :: pkk
    real(dp) :: v(2, 2)
    type(rng_state) :: state

    call rng_seed(state, 314159_8)
    call truncated_normal_sample(state, 0.0_dp, 1.0_dp, -1.0_dp, 0.5_dp, draw, info)
    if (info /= 0) error stop 'truncated-normal draw failed'

    v = reshape([1.0_dp, 0.2_dp, 0.2_dp, 1.5_dp], [2, 2])
    call riw_mcmcglmm(state, v, 6.0_dp, covariance, info)
    if (info /= 0) error stop 'inverse-Wishart draw failed'

    pkk = pkk_probability([0.5_dp, 0.5_dp], 2.0_dp)
    print '(a,f10.6)', 'truncated normal draw: ', draw
    print '(a,f10.6)', 'inverse-Wishart (1,1): ', covariance(1, 1)
    print '(a,f10.6)', 'pkk(0.5, 0.5, 2): ', pkk
end program distributions_example
