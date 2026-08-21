! SPDX-License-Identifier: GPL-2.0-only
program test_distributions
  use mvtnorm
  use test_support
  implicit none
  real(dp) :: mean(3),cov(3,3),x(3)
  real(dp),allocatable :: sim(:,:)
  integer :: j

  mean=[0.2_dp,-0.1_dp,0.3_dp]
  cov=reshape([1.0_dp,0.4_dp,-0.2_dp, 0.4_dp,1.5_dp,0.25_dp, -0.2_dp,0.25_dp,0.8_dp],[3,3])
  x=[0.5_dp,-0.4_dp,0.9_dp]
  call assert_close(dmvnorm_one(x,mean,cov,.true.),-3.249558631174606_dp,2.0e-13_dp,'dmvnorm log density')
  call assert_close(dmvt_one(x,mean,cov,7.0_dp,.true.),-3.3290726185114075_dp,2.0e-13_dp,'dmvt log density')
  call assert_close(normal_cdf(normal_quantile(0.12345_dp)),0.12345_dp,2.0e-14_dp,'normal inversion')
  call assert_close(student_t_cdf(student_t_quantile(0.91_dp,7.0_dp),7.0_dp),0.91_dp,2.0e-12_dp,'t inversion')
  call assert_close(chi_square_cdf(chi_square_quantile(0.73_dp,5.0_dp),5.0_dp),0.73_dp,2.0e-12_dp,'chi-square inversion')

  sim=rmvnorm(120000,mean,cov,24680)
  do j=1,3
    call assert_close(sum(sim(:,j))/real(size(sim,1),dp),mean(j),1.2e-2_dp,'rmvnorm mean')
  end do
  call assert_close(sum((sim(:,1)-mean(1))*(sim(:,2)-mean(2)))/real(size(sim,1)-1,dp),cov(1,2),1.8e-2_dp,'rmvnorm covariance')
  print '(a)', 'test_distributions: PASS'
end program test_distributions
