! SPDX-License-Identifier: GPL-3.0-only
program test_blmodel
  use blmodel, only : dp, posterior_result, observ_normal, bl_post_distribution
  use test_support, only : assert_true, assert_vector_close, assert_close
  implicit none

  real(dp) :: returns(4, 2), probabilities(4), weights(2), pick(2, 2), q(2), covariance(2, 2)
  type(posterior_result) :: result

  returns = reshape([0.01_dp, 0.03_dp, -0.02_dp, 0.00_dp, &
                     0.02_dp, -0.01_dp, 0.04_dp, 0.01_dp], [4, 2])
  probabilities = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
  weights = [0.6_dp, 0.4_dp]
  pick = 0.0_dp
  pick(1, 1) = 1.0_dp
  pick(2, 2) = 1.0_dp
  q = [0.024_dp, -0.012_dp]
  covariance = reshape([0.003708_dp, -0.003552_dp, -0.003552_dp, 0.003888_dp], [2, 2])

  result = bl_post_distribution(returns, probabilities, 12.0_dp, 'elliptic', weights, 0.5_dp, pick, q, &
    0.2_dp, 'MAD', 0.95_dp, observ_normal, 'diag', covariance)
  call assert_true(result%ok, 'Black-Litterman status')
  call assert_vector_close(result%equilibrium_returns, [0.00211030164096819_dp, -0.00151185789203691_dp], &
    3.0e-15_dp, 3.0e-14_dp, 'Black-Litterman equilibrium')
  call assert_vector_close(result%probabilities, [0.116321147284979_dp, 0.146780306959601_dp, &
    0.263508023626479_dp, 0.473390522128942_dp], 6.0e-14_dp, 6.0e-14_dp, &
    'Black-Litterman posterior')
  call assert_close(sum(result%probabilities), 1.0_dp, 2.0e-15_dp, 2.0e-14_dp, 'BL normalization')

  result = bl_post_distribution(returns, probabilities, 12.0_dp, 'bad-prior', weights, 0.5_dp, pick, q, &
    0.2_dp, 'MAD', 0.95_dp, observ_normal, 'diag', covariance)
  call assert_true(.not. result%ok, 'invalid prior rejected')

  write(*, '(a)') 'test_blmodel: PASS'
end program test_blmodel
