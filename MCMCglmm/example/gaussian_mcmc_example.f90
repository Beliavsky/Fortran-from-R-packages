! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program gaussian_mcmc_example
    use mcmcglmm, only : dp, rng_state, rng_seed, gaussian_mcmc_result, gaussian_mixed_mcmc
    implicit none

    integer :: i
    integer :: info
    real(dp) :: a_inverse(4, 4)
    real(dp) :: beta_prior_mean(1, 1)
    real(dp) :: beta_prior_precision(1, 1)
    real(dp) :: g_prior_scale(1, 1)
    real(dp) :: r_prior_scale(1, 1)
    real(dp) :: x(8, 1)
    real(dp) :: y(8, 1)
    real(dp) :: z(8, 4)
    type(gaussian_mcmc_result) :: result
    type(rng_state) :: state

    x(:, 1) = 1.0_dp
    y(:, 1) = [1.2_dp, 1.8_dp, 2.1_dp, 2.4_dp, 1.6_dp, 2.2_dp, 2.5_dp, 2.8_dp]
    z = 0.0_dp
    do i = 1, 8
        z(i, ceiling(real(i, dp) / 2.0_dp)) = 1.0_dp
    end do
    a_inverse = 0.0_dp
    do i = 1, 4
        a_inverse(i, i) = 1.0_dp
    end do
    beta_prior_mean = 0.0_dp
    beta_prior_precision = 0.01_dp
    g_prior_scale = 1.0_dp
    r_prior_scale = 1.0_dp

    call rng_seed(state, 12345_8)
    call gaussian_mixed_mcmc(y, x, z, a_inverse, beta_prior_mean, beta_prior_precision, &
        g_prior_scale, 4.0_dp, r_prior_scale, 4.0_dp, 600, 200, 4, state, result, info)
    if (info /= 0) error stop 'Gaussian mixed-model sampler failed'

    print '(a,f10.5)', 'posterior mean intercept: ', sum(result%beta(1, 1, :)) / real(size(result%beta, 3), dp)
    print '(a,f10.5)', 'posterior mean G variance: ', sum(result%g(1, 1, :)) / real(size(result%g, 3), dp)
end program gaussian_mcmc_example
