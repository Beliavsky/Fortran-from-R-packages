program family_engine_example
    use mcmcglmm, only : dp, rng_state, rng_seed, multivariate_family_mcmc_result, &
        ordinal_native_mcmc_result, heterogeneous_family_mixed_mcmc, ordinal_native_mixed_mcmc
    implicit none

    integer :: family(2)
    integer :: info
    integer :: ordinal_y(6)
    logical :: observed(6, 2)
    real(dp) :: a_inverse(2, 2)
    real(dp) :: additional(6, 2)
    real(dp) :: additional2(6, 2)
    real(dp) :: beta_mean(1, 2)
    real(dp) :: beta_mean_ordinal(1)
    real(dp) :: beta_precision(2, 2)
    real(dp) :: beta_precision_ordinal(1, 1)
    real(dp) :: g_scale(2, 2)
    real(dp) :: r_scale(2, 2)
    real(dp) :: x(6, 1)
    real(dp) :: y(6, 2)
    real(dp) :: z(6, 2)
    type(multivariate_family_mcmc_result) :: mixed_result
    type(ordinal_native_mcmc_result) :: ordinal_result
    type(rng_state) :: state

    x(:, 1) = 1.0_dp
    z = 0.0_dp
    z(1:3, 1) = 1.0_dp
    z(4:6, 2) = 1.0_dp
    a_inverse = 0.0_dp
    a_inverse(1, 1) = 1.0_dp
    a_inverse(2, 2) = 1.0_dp
    beta_mean = 0.0_dp
    beta_precision = 0.0_dp
    beta_precision(1, 1) = 0.1_dp
    beta_precision(2, 2) = 0.1_dp
    g_scale = 0.0_dp
    r_scale = 0.0_dp
    g_scale(1, 1) = 1.0_dp
    g_scale(2, 2) = 1.0_dp
    r_scale(1, 1) = 1.0_dp
    r_scale(2, 2) = 1.0_dp

    family = [1, 2]
    y(:, 1) = [-0.4_dp, 0.2_dp, 0.5_dp, -0.1_dp, 0.7_dp, 0.3_dp]
    y(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
    additional = 1.0_dp
    additional2 = 1.0_dp
    observed = .true.
    observed(2, 1) = .false.
    call rng_seed(state, 20260901_8)
    call heterogeneous_family_mixed_mcmc(family, y, additional, additional2, x, z, a_inverse, beta_mean, &
        beta_precision, g_scale, 5.0_dp, r_scale, 5.0_dp, .true., 0.55_dp, 120, 40, 4, state, mixed_result, info, &
        observed=observed)
    if (info /= 0) error stop 'heterogeneous sampler failed'
    print '(a,f8.4)', 'Gaussian-Poisson MH acceptance: ', mixed_result%acceptance_rate
    print '(a,2f10.5)', 'Last fixed-effect draw: ', mixed_result%beta(1, :, size(mixed_result%beta, 3))

    ordinal_y = [1, 1, 2, 2, 3, 3]
    beta_mean_ordinal = 0.0_dp
    beta_precision_ordinal(1, 1) = 0.1_dp
    call rng_seed(state, 20260902_8)
    call ordinal_native_mixed_mcmc(ordinal_y, [-1.0e40_dp, 0.0_dp, 1.0_dp, 1.0e40_dp], 0.25_dp, .true., &
        x, z, a_inverse, beta_mean_ordinal, beta_precision_ordinal, 1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, .true., &
        0.65_dp, 140, 40, 5, state, ordinal_result, info)
    if (info /= 0) error stop 'ordered-probit sampler failed'
    print '(a,f8.4)', 'Ordered-probit liability acceptance: ', ordinal_result%liability_acceptance_rate
    print '(a,f8.4)', 'Ordered-probit cutpoint acceptance: ', ordinal_result%cutpoint_acceptance_rate
end program family_engine_example
