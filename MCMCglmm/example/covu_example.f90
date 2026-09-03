program covu_example
    use mcmcglmm, only : dp, rng_state, rng_seed, covu_gaussian_mcmc_result, &
        covu_gaussian_mixed_mcmc
    implicit none

    integer :: info
    real(dp) :: beta_mean(1, 1)
    real(dp) :: beta_precision(1, 1)
    real(dp) :: initial_joint(2, 2)
    real(dp) :: joint_scale(2, 2)
    real(dp) :: loading(1, 1)
    real(dp) :: x(6, 1)
    real(dp) :: y(6, 1)
    type(covu_gaussian_mcmc_result) :: result
    type(rng_state) :: state

    x(:, 1) = 1.0_dp
    y(:, 1) = [-0.8_dp, -0.2_dp, 0.0_dp, 0.4_dp, 0.9_dp, 1.3_dp]
    loading(1, 1) = 1.0_dp
    beta_mean = 0.0_dp
    beta_precision(1, 1) = 0.2_dp
    joint_scale = reshape([1.0_dp, 0.15_dp, 0.15_dp, 1.0_dp], [2, 2])
    initial_joint = reshape([1.0_dp, 0.25_dp, 0.25_dp, 0.8_dp], [2, 2])

    call rng_seed(state, 817246_8)
    call covu_gaussian_mixed_mcmc(y, x, loading, beta_mean, beta_precision, &
        joint_scale, 4.0_dp, 48, 16, 4, state, result, info, initial_joint)
    if (info /= 0) error stop 'covu sampler failed'

    print '(a,f10.5)', 'mean beta: ', sum(result%beta(1, 1, :)) / real(size(result%beta, 3), dp)
    print '(a,f10.5)', 'mean marginal G: ', sum(result%g(1, 1, :)) / real(size(result%g, 3), dp)
    print '(a,f10.5)', 'mean conditional R: ', sum(result%r(1, 1, :)) / real(size(result%r, 3), dp)
    print '(a,f10.5)', 'mean residual-on-random regression: ', &
        sum(result%regression(1, 1, :)) / real(size(result%regression, 3), dp)
end program covu_example
