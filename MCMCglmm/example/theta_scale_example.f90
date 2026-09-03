program theta_scale_example
    use mcmcglmm, only : dp, rng_state, rng_seed, theta_scale_mcmc_result, &
        theta_scale_gaussian_mixed_mcmc
    implicit none

    integer :: info
    integer :: random_term(2)
    real(dp) :: a_inverse(2, 2)
    real(dp) :: beta_mean(1, 1)
    real(dp) :: beta_precision(1, 1)
    real(dp) :: g_df(1)
    real(dp) :: g_scale(1, 1, 1)
    real(dp) :: r_scale(1, 1)
    real(dp) :: x(6, 1)
    real(dp) :: x_scale(6, 1)
    real(dp) :: y(6, 1)
    real(dp) :: z(6, 2)
    real(dp) :: z_scale(6, 2)
    type(rng_state) :: state
    type(theta_scale_mcmc_result) :: result

    x(:, 1) = 1.0_dp
    x_scale(:, 1) = 1.0_dp
    z = 0.0_dp
    z(1:3, 1) = 1.0_dp
    z(4:6, 2) = 1.0_dp
    z_scale = 0.0_dp
    random_term = [1, 1]
    a_inverse = 0.0_dp
    a_inverse(1, 1) = 1.0_dp
    a_inverse(2, 2) = 1.0_dp
    beta_mean = 0.0_dp
    beta_precision(1, 1) = 0.2_dp
    g_scale(1, 1, 1) = 1.0_dp
    g_df(1) = 4.0_dp
    r_scale(1, 1) = 1.0_dp
    y(:, 1) = [-0.5_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.7_dp, 1.0_dp]

    call rng_seed(state, 913006_8)
    call theta_scale_gaussian_mixed_mcmc(y, x, z, x_scale, z_scale, random_term, &
        a_inverse, beta_mean, beta_precision, g_scale, g_df, r_scale, 4.0_dp, &
        1.0_dp, 1.0_dp, 48, 16, 4, state, result, info)
    if (info /= 0) error stop 'theta-scale sampler failed'

    print '(a,f10.5)', 'mean theta scale: ', sum(result%theta_scale) / real(size(result%theta_scale), dp)
    print '(a,f10.5)', 'mean G variance: ', sum(result%g(1, 1, 1, :)) / real(size(result%g, 4), dp)
    print '(a,f10.5)', 'mean R variance: ', sum(result%r(1, 1, :)) / real(size(result%r, 3), dp)
end program theta_scale_example
