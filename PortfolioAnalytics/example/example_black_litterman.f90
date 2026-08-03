! SPDX-License-Identifier: GPL-3.0-only
program example_black_litterman
  use portfolio_analytics
  implicit none
  real(dp) :: mu(3),sigma(3,3),pick(1,3),views(1),post_mu(3),post_sigma(3,3)
  integer :: info

  mu=[0.04_dp,0.05_dp,0.06_dp]
  sigma=reshape([0.04_dp,0.01_dp,0.00_dp,0.01_dp,0.05_dp,0.01_dp, &
                 0.00_dp,0.01_dp,0.08_dp],[3,3])
  pick(1,:)=[1.0_dp,-1.0_dp,0.0_dp]
  views=[0.03_dp]
  call black_litterman(mu,sigma,pick,views,post_mu,post_sigma,info=info)
  write(*,'(a,3f10.5)') 'prior means:',mu
  write(*,'(a,3f10.5)') 'posterior means:',post_mu
end program example_black_litterman
