! SPDX-License-Identifier: GPL-3.0-only
program test_posterior
  use blmodel, only : dp, posterior_result, observ_normal, post_distribution
  use test_support, only : assert_true, assert_vector_close, assert_matrix_close, assert_close
  implicit none

  real(dp) :: returns(4, 2), probabilities(4), equilibrium(2), q(2), pick(2, 2), covariance(2, 2)
  real(dp) :: expected_returns(4, 2), expected_covariance(2, 2)
  type(posterior_result) :: result

  returns = reshape([0.01_dp, 0.03_dp, -0.02_dp, 0.00_dp, &
                     0.02_dp, -0.01_dp, 0.04_dp, 0.01_dp], [4, 2])
  probabilities = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
  equilibrium = [0.005_dp, 0.01_dp]
  q = [0.002_dp, -0.001_dp]
  pick = 0.0_dp
  pick(1, 1) = 1.0_dp
  pick(2, 2) = 1.0_dp
  covariance = reshape([0.0004_dp, 0.0001_dp, 0.0001_dp, 0.0009_dp], [2, 2])

  result = post_distribution(returns, probabilities, equilibrium, q, pick, covariance, 0.2_dp, &
    observ_normal, 'diag')
  call assert_true(result%ok, 'posterior status')
  expected_returns = reshape([0.014_dp, 0.034_dp, -0.016_dp, 0.004_dp, &
                              0.014_dp, -0.016_dp, 0.034_dp, 0.004_dp], [4, 2])
  expected_covariance = 0.0_dp
  expected_covariance(1, 1) = 0.002_dp
  expected_covariance(2, 2) = 0.0045_dp
  call assert_matrix_close(result%returns, expected_returns, 2.0e-15_dp, 2.0e-14_dp, 'shifted returns')
  call assert_matrix_close(result%view_covariance, expected_covariance, 2.0e-15_dp, 2.0e-14_dp, &
    'view covariance')
  call assert_vector_close(result%probabilities, [0.106303879928644_dp, 0.170621723878166_dp, &
    0.272817469598847_dp, 0.450256926594344_dp], 5.0e-14_dp, 5.0e-14_dp, 'posterior probabilities')
  call assert_close(sum(result%probabilities), 1.0_dp, 2.0e-15_dp, 2.0e-14_dp, 'posterior normalization')

  write(*, '(a)') 'test_posterior: PASS'
end program test_posterior
