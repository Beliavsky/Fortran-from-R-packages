! SPDX-License-Identifier: GPL-2.0-only
program clustering_mmd
  use kernlab
  implicit none
  real(dp) :: x(8,2), y(8,2)
  type(kernel_spec) :: kernel
  type(cluster_result) :: clusters
  type(mmd_result) :: test

  x = reshape([ -2.0_dp,-2.0_dp, -1.8_dp,-2.2_dp, -2.2_dp,-1.8_dp, -1.9_dp,-1.7_dp, &
                 2.0_dp, 2.0_dp,  1.8_dp, 2.2_dp,  2.2_dp, 1.8_dp,  1.9_dp, 1.7_dp ], &
               [8,2], order=[2,1])
  y = x + 0.75_dp
  kernel = rbfdot(0.5_dp)

  call specc(x, 2, kernel, clusters)
  print '(a,8(1x,i0))', 'Spectral-cluster labels:', clusters%labels

  call kmmd(x, y, kernel, test, bootstrap=.true., ntimes=50)
  print '(a,f10.6)', 'MMD statistic: ', test%mmd1
  print '(a,l1)', 'Reject by Rademacher bound: ', test%reject_rademacher
end program clustering_mmd
