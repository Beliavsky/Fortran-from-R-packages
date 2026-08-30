! SPDX-License-Identifier: GPL-2.0-only
program test_probabilities
  use mvtnorm
  use test_support
  implicit none
  real(dp) :: mean(3),cov(3,3),lo(3),up(3),cor2(2,2),lo2(2),up2(2)
  type(probability_control) :: ctl
  type(probability_result) :: pr

  call assert_close(bivariate_normal_cdf(0.3_dp,0.7_dp,0.5_dp),0.5328056535231198_dp,2.0e-12_dp,'bivariate cdf')
  cor2=reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
  lo2=[-huge(1.0_dp),-huge(1.0_dp)]; up2=[0.3_dp,0.7_dp]
  ctl%method=method_tvpack
  pr=rectangle_probability_normal(lo2,up2,cor2,ctl)
  call assert_close(pr%value,0.5328056535231198_dp,2.0e-12_dp,'TVPACK bivariate route')

  mean=[0.2_dp,-0.1_dp,0.3_dp]
  cov=reshape([1.0_dp,0.4_dp,-0.2_dp, 0.4_dp,1.5_dp,0.25_dp, -0.2_dp,0.25_dp,0.8_dp],[3,3])
  lo=[-1.0_dp,-0.5_dp,-0.8_dp]; up=[0.7_dp,1.2_dp,0.6_dp]
  ctl=probability_control(maxpts=400000,batches=20,abseps=2.0e-5_dp,releps=0.0_dp,seed=12345)
  pr=pmvnorm(lo,up,mean,cov,ctl)
  call assert_close(pr%value,0.1388616650_dp,2.0e-5_dp,'trivariate normal rectangle')
  call assert_true(pr%error<5.0e-5_dp,'normal error estimate')
  pr=pmvt(lo,up,mean,cov,7.0_dp,ctl)
  call assert_close(pr%value,0.13325511_dp,7.0e-5_dp,'trivariate t rectangle')

  ctl%method=method_miwa; ctl%maxpts=120000; ctl%batches=20
  pr=pmvnorm(lo,up,mean,cov,ctl)
  call assert_close(pr%value,0.1388616650_dp,8.0e-5_dp,'Miwa-compatible route')
  print '(a)', 'test_probabilities: PASS'
end program test_probabilities
