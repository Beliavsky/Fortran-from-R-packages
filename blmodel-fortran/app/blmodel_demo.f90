! SPDX-License-Identifier: GPL-3.0-only
program blmodel_demo
  use blmodel, only : dp, posterior_result, observ_normal, bl_post_distribution
  implicit none

  real(dp) :: returns(4, 2), probabilities(4), weights(2), pick(2, 2), q(2)
  type(posterior_result) :: result
  integer :: i

  returns = reshape([0.01_dp, 0.03_dp, -0.02_dp, 0.00_dp, &
                     0.02_dp, -0.01_dp, 0.04_dp, 0.01_dp], [4, 2])
  probabilities = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp]
  weights = [0.6_dp, 0.4_dp]
  pick = 0.0_dp
  pick(1, 1) = 1.0_dp
  pick(2, 2) = 1.0_dp
  q = [0.024_dp, -0.012_dp]

  result = bl_post_distribution(returns, probabilities, 12.0_dp, 'elliptic', weights, 0.5_dp, pick, q, &
    0.2_dp, 'MAD', 0.95_dp, observ_normal, 'diag')
  if (.not. result%ok) error stop result%message

  write(*, '(a,2(1x,f12.8))') 'equilibrium returns:', result%equilibrium_returns
  write(*, '(a)') 'posterior scenarios:'
  do i = 1, size(result%probabilities)
    write(*, '(i3,2(1x,f10.6),1x,f10.6)') i, result%returns(i, :), result%probabilities(i)
  end do
end program blmodel_demo
