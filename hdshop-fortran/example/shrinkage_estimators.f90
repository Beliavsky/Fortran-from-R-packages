! SPDX-License-Identifier: GPL-3.0-only
program shrinkage_estimators
  use hdshop, only: dp, sigma_sample_estimator, cov_shrink_bgp14, &
    matrix_shrink_result, mean_js, mean_shrink_result
  implicit none
  real(dp)::x(2,6),target(2,2)
  real(dp),allocatable::s(:,:)
  type(matrix_shrink_result)::covfit
  type(mean_shrink_result)::meanfit
  x=reshape([0.01_dp,0.02_dp,0.00_dp,0.01_dp,-0.01_dp,0.03_dp, &
    0.02_dp,0.00_dp,0.01_dp,-0.02_dp,0.03_dp,0.01_dp],[2,6])
  s=sigma_sample_estimator(x);target=0.0_dp;target(1,1)=0.5_dp;target(2,2)=0.5_dp
  covfit=cov_shrink_bgp14(6,target,s);meanfit=mean_js(x,0.0_dp)
  print '(a,f10.6)', 'covariance alpha: ',covfit%alpha
  print '(a,2f10.6)', 'shrunk means: ',meanfit%means
end program shrinkage_estimators
