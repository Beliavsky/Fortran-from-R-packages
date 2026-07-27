! SPDX-License-Identifier: GPL-3.0-only
program test_equilibrium
  use blmodel, only : dp, moment_result, equilibrium_result, discrete_variance, equilibrium_mean, &
    equilibrium_mean_elliptic, diag_of, make_diag
  use test_support, only : assert_true, assert_vector_close, assert_matrix_close, assert_close
  implicit none

  real(dp) :: returns(4, 2), probabilities(4), weights(2), covariance(2, 2)
  type(moment_result) :: moments
  type(equilibrium_result) :: equilibrium

  returns = reshape([0.01_dp, 0.03_dp, -0.02_dp, 0.00_dp, &
                     0.02_dp, -0.01_dp, 0.04_dp, 0.01_dp], [4, 2])
  probabilities = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
  weights = [0.6_dp, 0.4_dp]

  moments = discrete_variance(returns, probabilities, 12.0_dp)
  call assert_true(moments%ok, 'discrete moments status')
  call assert_vector_close(moments%mean, [0.012_dp, 0.192_dp], 1.0e-15_dp, 1.0e-14_dp, 'annual mean')
  covariance = reshape([0.003708_dp, -0.003552_dp, -0.003552_dp, 0.003888_dp], [2, 2])
  call assert_matrix_close(moments%covariance, covariance, 2.0e-15_dp, 2.0e-14_dp, 'annual covariance')
  call assert_vector_close(diag_of(covariance), [0.003708_dp, 0.003888_dp], 0.0_dp, 0.0_dp, 'diag_of')
  call assert_matrix_close(make_diag([2.0_dp, 3.0_dp]), reshape([2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp], [2, 2]), &
    0.0_dp, 0.0_dp, 'make_diag')

  equilibrium = equilibrium_mean(returns, probabilities, weights, 0.01_dp, 'MAD')
  call assert_true(equilibrium%ok, 'MAD equilibrium status')
  call assert_vector_close(equilibrium%market_returns, [0.0319047619047619_dp, -0.0228571428571429_dp], &
    2.0e-14_dp, 2.0e-14_dp, 'MAD equilibrium')
  call assert_close(dot_product(equilibrium%market_returns, weights), 0.01_dp, 2.0e-15_dp, 2.0e-14_dp, &
    'MAD market return')

  equilibrium = equilibrium_mean(returns, probabilities, weights, 0.01_dp, 'CVAR', 0.75_dp)
  call assert_true(equilibrium%ok, 'CVAR equilibrium status')
  call assert_vector_close(equilibrium%market_returns, [0.00333333333333333_dp, 0.02_dp], &
    2.0e-14_dp, 2.0e-14_dp, 'CVAR equilibrium')

  equilibrium = equilibrium_mean_elliptic(covariance / 12.0_dp, weights, &
    0.5_dp / sqrt(dot_product(weights, matmul(covariance, weights))))
  call assert_true(equilibrium%ok, 'elliptic equilibrium status')
  call assert_vector_close(equilibrium%market_returns, [0.00211030164096819_dp, -0.00151185789203691_dp], &
    3.0e-15_dp, 3.0e-14_dp, 'elliptic equilibrium')

  write(*, '(a)') 'test_equilibrium: PASS'
end program test_equilibrium
