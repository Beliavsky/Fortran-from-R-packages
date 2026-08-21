! SPDX-License-Identifier: GPL-2.0-only
program distributions_and_probabilities
  use mvtnorm
  implicit none
  real(dp) :: mean(2),sigma(2,2),lower(2),upper(2)
  real(dp),allocatable :: draws(:,:)
  type(probability_result) :: pn,pt
  type(probability_control) :: control

  mean=[0.0_dp,0.0_dp]
  sigma=reshape([1.0_dp,0.6_dp,0.6_dp,1.0_dp],[2,2])
  lower=[-1.0_dp,-1.0_dp]
  upper=[ 1.0_dp, 1.0_dp]
  control=tvpack(1.0e-8_dp)

  pn=pmvnorm(lower,upper,mean,sigma,control)
  control=genz_bretz(maxpts=160000,abseps=2.0e-4_dp,seed=42)
  pt=pmvt(lower,upper,mean,sigma,5.0_dp,control)
  draws=rmvt_shifted(5,sigma,5.0_dp,mean,42)

  print '(a,f12.8)', 'normal probability: ',pn%value
  print '(a,f12.8)', 't probability:      ',pt%value
  print '(a)', 'five t draws:'
  print '(2f12.6)',transpose(draws)
end program distributions_and_probabilities
