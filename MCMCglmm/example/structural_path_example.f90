program structural_path_example
    use mcmcglmm, only : dp, rng_state, rng_seed, structural_gaussian_mcmc_result, &
        structural_gaussian_multi_term_mcmc
    implicit none

    integer :: i
    integer :: info
    integer :: random_term(2)
    real(dp) :: a_inverse(2, 2)
    real(dp) :: basis(2, 2, 1)
    real(dp) :: beta_mean(1, 2)
    real(dp) :: beta_precision(2, 2)
    real(dp) :: g_df(1)
    real(dp) :: g_scale(2, 2, 1)
    real(dp) :: r_scale(2, 2)
    real(dp) :: structural_mean(1)
    real(dp) :: structural_precision(1, 1)
    real(dp) :: structural_sd(1)
    real(dp) :: x(8, 1)
    real(dp) :: y(8, 2)
    real(dp) :: z(8, 2)
    type(rng_state) :: state
    type(structural_gaussian_mcmc_result) :: result

    x(:, 1) = 1.0_dp
    z = 0.0_dp
    z(1:4, 1) = 1.0_dp
    z(5:8, 2) = 1.0_dp
    random_term = [1, 1]
    a_inverse = 0.0_dp
    do i = 1, 2
        a_inverse(i, i) = 1.0_dp
    end do
    y(:, 1) = [-0.8_dp, -0.4_dp, -0.1_dp, 0.2_dp, 0.5_dp, 0.7_dp, 1.0_dp, 1.3_dp]
    y(:, 2) = [-0.2_dp, -0.1_dp, 0.1_dp, 0.3_dp, 0.5_dp, 0.8_dp, 1.0_dp, 1.2_dp]
    beta_mean = 0.0_dp
    beta_precision = 0.0_dp
    beta_precision(1, 1) = 0.2_dp
    beta_precision(2, 2) = 0.2_dp
    g_scale = 0.0_dp
    g_scale(1, 1, 1) = 1.0_dp
    g_scale(2, 2, 1) = 1.0_dp
    g_df(1) = 5.0_dp
    r_scale = 0.0_dp
    r_scale(1, 1) = 1.0_dp
    r_scale(2, 2) = 1.0_dp
    basis = 0.0_dp
    basis(2, 1, 1) = 1.0_dp
    structural_mean = 0.0_dp
    structural_precision(1, 1) = 0.5_dp
    structural_sd(1) = 0.08_dp

    call rng_seed(state, 913007_8)
    call structural_gaussian_multi_term_mcmc(y, x, z, random_term, a_inverse, &
        beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, basis, &
        structural_mean, structural_precision, structural_sd, 60, 20, 4, state, &
        result, info)
    if (info /= 0) error stop 'structural path sampler failed'

    print '(a,f10.5)', 'mean path coefficient: ', &
        sum(result%structural_parameter(1, :)) / real(size(result%structural_parameter, 2), dp)
    print '(a,f10.5)', 'path MH acceptance: ', result%structural_acceptance_rate
end program structural_path_example
