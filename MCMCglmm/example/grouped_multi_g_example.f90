program grouped_multi_g_example
    use mcmcglmm, only : dp, rng_state, rng_seed, grouped_multi_term_mcmc_result, &
        two_part_multi_term_mixed_mcmc
    implicit none

    integer :: i
    integer :: info
    integer :: random_term(4)
    integer :: trials(8)
    logical :: observed(8)
    real(dp) :: a_inverse(4, 4)
    real(dp) :: beta_mean(1, 2)
    real(dp) :: beta_precision(2, 2)
    real(dp) :: g_df(2)
    real(dp) :: g_scale(2, 2, 2)
    real(dp) :: r_scale(2, 2)
    real(dp) :: x(8, 1)
    real(dp) :: y(8)
    real(dp) :: z(8, 4)
    type(grouped_multi_term_mcmc_result) :: result
    type(rng_state) :: state

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
    y = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp, 4.0_dp]
    trials = 1
    observed = .true.
    observed(3) = .false.

    call rng_seed(state, 817236_8)
    call two_part_multi_term_mixed_mcmc(11, y, trials, x, z, random_term, &
        a_inverse, beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, &
        .true., 0.30_dp, 48, 16, 4, state, result, info, observed=observed)
    if (info /= 0) error stop 'grouped multi-G sampler failed'

    print '(a,f10.5)', 'latent MH acceptance: ', result%acceptance_rate
    print '(a,2f10.5)', 'mean first-trait G variances: ', &
        sum(result%g(1, 1, 1, :)) / real(size(result%g, 4), dp), &
        sum(result%g(1, 1, 2, :)) / real(size(result%g, 4), dp)
    print '(a,l1)', 'missing row liability finite: ', all(abs(result%last_liability(3, :)) < huge(1.0_dp))
end program grouped_multi_g_example
