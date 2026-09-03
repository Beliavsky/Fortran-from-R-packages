! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program unified_multi_g_example
    use mcmcglmm, only : dp, heterogeneous_multi_term_mixed_mcmc, rng_seed, rng_state, unified_family_mcmc_result
    implicit none

    integer :: family(2)
    integer :: i
    integer :: info
    integer :: random_term(4)
    logical :: observed(8, 2)
    real(dp) :: additional(8, 2)
    real(dp) :: additional2(8, 2)
    real(dp) :: a_inverse(4, 4)
    real(dp) :: beta_mean(1, 2)
    real(dp) :: beta_precision(2, 2)
    real(dp) :: g_df(2)
    real(dp) :: g_scale(2, 2, 2)
    real(dp) :: r_scale(2, 2)
    real(dp) :: x(8, 1)
    real(dp) :: y(8, 2)
    real(dp) :: z(8, 4)
    type(unified_family_mcmc_result) :: result
    type(rng_state) :: state

    family = [1, 2]
    x(:, 1) = 1.0_dp
    z = 0.0_dp
    z(1:4, 1) = 1.0_dp
    z(5:8, 2) = 1.0_dp
    z([1, 3, 5, 7], 3) = 1.0_dp
    z([2, 4, 6, 8], 4) = 1.0_dp
    random_term = [1, 1, 2, 2]
    a_inverse = 0.0_dp
    do i = 1, 4
        a_inverse(i, i) = 1.0_dp
    end do
    y(:, 1) = [-0.4_dp, -0.1_dp, 0.1_dp, 0.2_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.1_dp]
    y(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp]
    observed = .true.
    observed(2, 2) = .false.
    y(2, 2) = -999.0_dp
    additional = 0.0_dp
    additional2 = 0.0_dp
    beta_mean = 0.0_dp
    beta_precision = 0.0_dp
    beta_precision(1, 1) = 0.1_dp
    beta_precision(2, 2) = 0.1_dp
    g_scale = 0.0_dp
    g_scale(1, 1, :) = 1.0_dp
    g_scale(2, 2, :) = 1.0_dp
    g_df = 5.0_dp
    r_scale = 0.0_dp
    r_scale(1, 1) = 1.0_dp
    r_scale(2, 2) = 1.0_dp

    call rng_seed(state, 20260904_8)
    call heterogeneous_multi_term_mixed_mcmc(family, y, additional, additional2, x, z, random_term, a_inverse, &
        beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.35_dp, &
        100, 40, 4, state, result, info, observed)
    if (info /= 0) error stop 'unified mixed-family multi-G sampler failed'

    print '(a,f8.4)', 'latent MH acceptance: ', result%acceptance_rate
    print '(a,2f10.5)', 'last fixed-effect draw: ', result%beta(1, :, size(result%beta, 3))
end program unified_multi_g_example
