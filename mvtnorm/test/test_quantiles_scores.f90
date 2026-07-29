! SPDX-License-Identifier: GPL-2.0-only
program test_quantiles_scores
  use mvtnorm
  use test_support
  implicit none
  real(dp) :: mean(2),cov(2,2),obs(1,2),lower(2),upper(2),b(2,1),d(2),z(1,3),w(3)
  real(dp) :: expected_q,rr
  type(probability_control) :: ctl
  type(quantile_result) :: qr
  type(likelihood_result) :: sr

  mean=0.0_dp; cov=0.0_dp; cov(1,1)=1.0_dp; cov(2,2)=1.0_dp
  ctl=probability_control(maxpts=60000,batches=12,abseps=1.0e-5_dp,seed=123)
  qr=qmvnorm(0.81_dp,mean,cov,'lower',ctl,lower_bound=-2.0_dp,upper_bound=3.0_dp,ptol=2.0e-5_dp)
  expected_q=normal_quantile(sqrt(0.81_dp))
  call assert_true(qr%converged,'qmvnorm convergence')
  call assert_close(qr%quantile,expected_q,4.0e-5_dp,'qmvnorm independent reference')

  obs(1,:)=[0.4_dp,-0.7_dp]
  sr=sldmvnorm(obs,mean,cov,1.0e-5_dp)
  call assert_true(sr%ok,'density score')
  call assert_close(sr%score(1,1),0.4_dp,2.0e-7_dp,'mean score 1')
  call assert_close(sr%score(2,1),-0.7_dp,2.0e-7_dp,'mean score 2')

  lower=[-1.0_dp,-0.5_dp]; upper=[0.8_dp,1.1_dp]; mean=[0.1_dp,-0.2_dp]
  b(:,1)=[0.4_dp,-0.3_dp]; d=[0.7_dp,1.2_dp]
  z=reshape([-1.0_dp,0.0_dp,1.0_dp],[1,3]); w=1.0_dp/3.0_dp
  rr=lpRR(lower,upper,mean,b,d,z,w,.false.)
  call assert_close(rr,0.3250046942035033_dp,2.0e-6_dp,'random-effects probability')
  print '(a)', 'test_quantiles_scores: PASS'
end program test_quantiles_scores
