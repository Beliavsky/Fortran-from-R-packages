! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program threshold_mcmc_example
    use mcmcglmm, only : dp, rng_state, rng_seed, ordinal_mcmc_result, ordinal_probit_mixed_mcmc
    implicit none

    integer :: i
    integer :: info
    integer :: y(8)
    real(dp) :: a_inverse(4, 4)
    real(dp) :: beta_mean(1)
    real(dp) :: beta_precision(1, 1)
    real(dp) :: cutpoints(3)
    real(dp) :: x(8, 1)
    real(dp) :: z(8, 4)
    type(ordinal_mcmc_result) :: result
    type(rng_state) :: state

    y = [1, 1, 1, 2, 1, 2, 2, 2]
    cutpoints = [-1.0e40_dp, 0.0_dp, 1.0e40_dp]
    x(:, 1) = 1.0_dp
    z = 0.0_dp
    do i = 1, 8
        z(i, ceiling(real(i, dp) / 2.0_dp)) = 1.0_dp
    end do
    a_inverse = 0.0_dp
    do i = 1, 4
        a_inverse(i, i) = 1.0_dp
    end do
    beta_mean = 0.0_dp
    beta_precision = 0.01_dp

    call rng_seed(state, 987654_8)
    call ordinal_probit_mixed_mcmc(y, cutpoints, x, z, a_inverse, beta_mean, beta_precision, &
        1.0_dp, 4.0_dp, 600, 200, 4, state, result, info)
    if (info /= 0) error stop 'threshold-probit sampler failed'

    print '(a,f10.5)', 'posterior mean intercept: ', sum(result%beta(1, :)) / real(size(result%beta, 2), dp)
    print '(a,f10.5)', 'posterior mean random variance: ', sum(result%g) / real(size(result%g), dp)
end program threshold_mcmc_example
