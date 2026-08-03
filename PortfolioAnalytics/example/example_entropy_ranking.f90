! SPDX-License-Identifier: GPL-3.0-only
program example_entropy_ranking
  use portfolio_analytics
  implicit none
  real(dp) :: returns(6,3),prior(6),mu(3),sigma(3,3),posterior(6)
  integer :: order(3)
  logical :: converged

  returns(:,1)=[-0.03_dp,-0.01_dp,0.00_dp,0.01_dp,0.02_dp,0.04_dp]
  returns(:,2)=[0.02_dp,0.01_dp,0.00_dp,-0.01_dp,-0.02_dp,-0.03_dp]
  returns(:,3)=[-0.01_dp,0.00_dp,0.01_dp,0.02_dp,0.03_dp,0.04_dp]
  prior=1.0_dp/6.0_dp
  order=[2,1,3]
  call meucci_ranking(returns,prior,order,mu,sigma,posterior,converged)
  write(*,'(a,l1)') 'converged:',converged
  write(*,'(a,3f10.5)') 'rank-constrained means:',mu
end program example_entropy_ranking
