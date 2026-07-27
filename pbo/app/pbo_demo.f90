! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
program pbo_demo
  use pbo, only : dp, pbo_result, compute_pbo, sharpe_ratio
  implicit none
  type(pbo_result) :: result
  real(dp) :: returns(12,4)
  integer :: i

  do i = 1, 12
    returns(i,1) = 0.001_dp * real(i,dp)
    returns(i,2) = 0.008_dp - 0.0008_dp * real(i,dp)
    returns(i,3) = 0.003_dp * sin(real(i,dp))
    returns(i,4) = merge(0.006_dp, -0.004_dp, mod(i,2) == 0)
  end do
  call compute_pbo(returns, 4, sharpe_ratio, result)
  if (.not. result%success) error stop result%message
  print '(a,f8.4)', 'PBO: ', result%phi
  print '(a,f8.4)', 'Probability of loss: ', result%below_threshold
  print '(a,f10.5)', 'Degradation slope: ', result%degradation_slope
end program pbo_demo
