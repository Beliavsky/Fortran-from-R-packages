! SPDX-License-Identifier: GPL-2.0-only
program test_conditioning_likelihood
  use mvtnorm
  use test_support
  implicit none
  real(dp) :: mean(3),cov(3,3),given(2),obs(2,3),low(2,3),upp(2,3)
  real(dp) :: exact(2,1),mlow(2,2),mupp(2,2)
  integer :: idx(2)
  type(conditional_result) :: cr
  type(likelihood_result) :: lr
  type(probability_control) :: ctl

  mean=[0.2_dp,-0.1_dp,0.3_dp]
  cov=reshape([1.0_dp,0.4_dp,-0.2_dp, 0.4_dp,1.5_dp,0.25_dp, -0.2_dp,0.25_dp,0.8_dp],[3,3])
  idx=[1,3]; given=[0.4_dp,0.1_dp]
  cr=conditional_mvnormal(mean,cov,idx,given)
  call assert_true(cr%ok,'conditional result')
  call assert_close(cr%mean(1),-0.0894736842105263_dp,2.0e-13_dp,'conditional mean')
  call assert_close(cr%covariance(1,1),1.1967105263157894_dp,2.0e-13_dp,'conditional covariance')

  obs=reshape([0.5_dp,0.1_dp,-0.4_dp,0.2_dp,0.9_dp,-0.2_dp],[2,3])
  lr=ldmvnorm(obs,mean,cov,.true.)
  call assert_close(lr%loglik(1),-3.249558631174606_dp,2.0e-13_dp,'exact likelihood 1')
  call assert_close(lr%loglik(2),-3.0260292194098994_dp,2.0e-10_dp,'exact likelihood 2')

  low=reshape([-1.0_dp,-0.2_dp,-0.5_dp,-0.4_dp,-0.8_dp,-0.1_dp],[2,3])
  upp=reshape([0.7_dp,1.0_dp,1.2_dp,0.8_dp,0.6_dp,1.2_dp],[2,3])
  ctl=probability_control(maxpts=220000,batches=16,abseps=3.0e-5_dp,seed=9981)
  lr=lpmvnorm(low,upp,mean,cov,ctl,.true.)
  call assert_close(exp(lr%loglik(1)),0.1388617_dp,7.0e-5_dp,'interval likelihood 1')
  call assert_close(exp(lr%loglik(2)),0.09362_dp,7.0e-5_dp,'interval likelihood 2')

  exact(:,1)=[0.3_dp,-0.2_dp]
  mlow=reshape([-0.5_dp,-0.2_dp,-0.6_dp,-0.3_dp],[2,2])
  mupp=reshape([0.8_dp,1.0_dp,0.7_dp,0.9_dp],[2,2])
  lr=ldpmvnorm(exact,mlow,mupp,mean,cov,ctl,.true.)
  call assert_close(lr%loglik(1),-2.4054622_dp,2.0e-4_dp,'mixed likelihood 1')
  call assert_close(lr%loglik(2),-2.7387028_dp,2.0e-4_dp,'mixed likelihood 2')
  print '(a)', 'test_conditioning_likelihood: PASS'
end program test_conditioning_likelihood
