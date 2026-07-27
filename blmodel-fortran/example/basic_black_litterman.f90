! SPDX-License-Identifier: GPL-3.0-only
program basic_black_litterman
  use blmodel, only : dp, posterior_result, observ_student_t, bl_post_distribution
  implicit none

  real(dp) :: returns(5, 3), probabilities(5), weights(3), pick(2, 3), q(2), params(1)
  type(posterior_result) :: result

  returns = reshape([0.02_dp, -0.01_dp, 0.03_dp, 0.01_dp, -0.02_dp, &
                     0.01_dp, 0.04_dp, -0.02_dp, 0.00_dp, 0.03_dp, &
                    -0.01_dp, 0.02_dp, 0.01_dp, 0.04_dp, 0.00_dp], [5, 3])
  probabilities = 0.2_dp
  weights = [0.4_dp, 0.35_dp, 0.25_dp]
  pick = reshape([1.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp], [2, 3])
  q = [0.03_dp, 0.01_dp]
  params = 7.0_dp

  result = bl_post_distribution(returns, probabilities, 12.0_dp, 'elliptic', weights, 0.4_dp, pick, q, &
    0.05_dp, 'MAD', 0.95_dp, observ_student_t, 'full', view_params=params)
  if (.not. result%ok) error stop result%message
  write(*, '(a,3(1x,f11.7))') 'equilibrium:', result%equilibrium_returns
  write(*, '(a,5(1x,f10.6))') 'posterior probabilities:', result%probabilities
end program basic_black_litterman
