! SPDX-License-Identifier: GPL-2.0-only
program mvtnorm_demo
  use mvtnorm
  implicit none
  real(dp) :: mean(3),sigma(3,3),lower(3),upper(3),x(3)
  type(probability_control) :: control
  type(probability_result) :: probability
  type(quantile_result) :: quantile

  mean=[0.0_dp,0.2_dp,-0.1_dp]
  sigma=reshape([1.0_dp,0.5_dp,0.2_dp, 0.5_dp,1.5_dp,-0.1_dp, 0.2_dp,-0.1_dp,0.8_dp],[3,3])
  lower=[-1.0_dp,-0.8_dp,-0.5_dp]
  upper=[ 1.0_dp, 1.2_dp, 0.9_dp]
  x=[0.1_dp,0.4_dp,-0.2_dp]
  control=genz_bretz(maxpts=120000,abseps=2.0e-4_dp,seed=2026)

  probability=pmvnorm(lower,upper,mean,sigma,control)
  quantile=qmvnorm(0.90_dp,mean,sigma,'lower',control,ptol=5.0e-4_dp)

  print '(a,f12.8)', 'log density: ',dmvnorm_one(x,mean,sigma,.true.)
  print '(a,f12.8,a,es10.2)', 'rectangle probability: ',probability%value,' +/- ',probability%error
  print '(a,f12.8)', 'joint lower-tail 90% quantile: ',quantile%quantile
end program mvtnorm_demo
