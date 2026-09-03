! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program parameter_expansion_example
    use mcmcglmm, only : dp, gaussian_px_mcmc_result, gaussian_parameter_expanded_mcmc, rng_seed, rng_state
    implicit none

    integer :: info
    real(dp) :: a_inverse(2, 2)
    real(dp) :: alpha_precision(1, 1)
    real(dp) :: beta_mean(1, 1)
    real(dp) :: beta_precision(1, 1)
    real(dp) :: g_scale(1, 1)
    real(dp) :: r_scale(1, 1)
    real(dp) :: x(4, 1)
    real(dp) :: y(4, 1)
    real(dp) :: z(4, 2)
    type(gaussian_px_mcmc_result) :: result
    type(rng_state) :: state

    x(:, 1) = 1.0_dp
    y(:, 1) = [-0.4_dp, 0.2_dp, 0.8_dp, 0.5_dp]
    z = 0.0_dp
    z(1:2, 1) = 1.0_dp
    z(3:4, 2) = 1.0_dp
    a_inverse = 0.0_dp
    a_inverse(1, 1) = 1.0_dp
    a_inverse(2, 2) = 1.0_dp
    beta_mean = 0.0_dp
    beta_precision(1, 1) = 0.1_dp
    g_scale(1, 1) = 1.0_dp
    r_scale(1, 1) = 1.0_dp
    alpha_precision(1, 1) = 1.0_dp

    call rng_seed(state, 20260903_8)
    call gaussian_parameter_expanded_mcmc(y, x, z, a_inverse, beta_mean, beta_precision, g_scale, 4.0_dp, &
        r_scale, 4.0_dp, [1.0_dp], alpha_precision, .true., 120, 40, 4, state, result, info)
    if (info /= 0) error stop 'parameter-expanded Gaussian sampler failed'

    print '(a,f10.5)', 'posterior mean alpha: ', sum(result%alpha(1, :)) / real(size(result%alpha, 2), dp)
    print '(a,f10.5)', 'posterior mean expanded G: ', sum(result%g(1, 1, :)) / real(size(result%g, 3), dp)
end program parameter_expansion_example
