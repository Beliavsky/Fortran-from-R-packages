program prediction_example
    use mcmcglmm, only : dp, multi_term_build_v, posterior_linear_predictor, &
        scalar_response_expectation
    implicit none

    integer :: info
    integer :: random_term(2)
    logical :: include_g(2)
    logical :: marginalize_g(2)
    real(dp) :: a_inverse(2, 2)
    real(dp) :: beta_draws(1, 2, 1)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: diagonal(:, :)
    real(dp) :: expectation
    real(dp) :: g_matrix(2, 2, 2)
    real(dp), allocatable :: predictor(:, :, :)
    real(dp) :: r_matrix(2, 2)
    real(dp) :: random_draws(2, 2, 1)
    real(dp) :: x(2, 1)
    real(dp) :: z(2, 2)

    z = 0.0_dp
    z(1, 1) = 1.0_dp
    z(2, 2) = 1.0_dp
    random_term = [1, 2]
    a_inverse = 0.0_dp
    a_inverse(1, 1) = 1.0_dp
    a_inverse(2, 2) = 1.0_dp
    g_matrix = 0.0_dp
    g_matrix(:, :, 1) = reshape([1.0_dp, 0.2_dp, 0.2_dp, 2.0_dp], [2, 2])
    g_matrix(:, :, 2) = reshape([0.5_dp, 0.1_dp, 0.1_dp, 0.7_dp], [2, 2])
    r_matrix = reshape([0.3_dp, 0.05_dp, 0.05_dp, 0.4_dp], [2, 2])
    include_g = .true.

    call multi_term_build_v(z, random_term, a_inverse, g_matrix, r_matrix, &
        include_g, covariance, diagonal, info)
    if (info /= 0) error stop 'buildV failed'
    print '(a,2f10.4)', 'marginal variances at row 1: ', diagonal(1, :)

    x(:, 1) = 1.0_dp
    beta_draws(1, :, 1) = [1.0_dp, 10.0_dp]
    random_draws(:, :, 1) = reshape([0.5_dp, 2.0_dp, 1.5_dp, 3.0_dp], [2, 2])
    marginalize_g = [.false., .true.]
    call posterior_linear_predictor(x, z, random_term, beta_draws, random_draws, &
        marginalize_g, predictor, info)
    if (info /= 0) error stop 'posterior predictor failed'
    print '(a,2f10.4)', 'row 1 linear predictor: ', predictor(1, :, 1)

    call scalar_response_expectation(2, 0.3_dp, 0.4_dp, 0.0_dp, 0.0_dp, expectation, info)
    if (info /= 0) error stop 'response expectation failed'
    print '(a,f10.5)', 'Poisson response mean: ', expectation
end program prediction_example
