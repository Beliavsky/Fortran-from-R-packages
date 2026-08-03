! SPDX-License-Identifier: GPL-2.0-only
program demo_kernlab
  use kernlab
  implicit none
  real(dp) :: x(8,2)
  integer :: y(8), status
  type(kernel_spec) :: kernel
  type(kernel_model) :: model
  type(cluster_result) :: clusters
  real(dp), allocatable :: scores(:,:)
  integer, allocatable :: predicted(:)

  x = reshape([ -2.0_dp,-2.0_dp, -1.8_dp,-2.2_dp, -2.2_dp,-1.8_dp, -1.9_dp,-1.7_dp, &
                 2.0_dp, 2.0_dp,  1.8_dp, 2.2_dp,  2.2_dp, 1.8_dp,  1.9_dp, 1.7_dp ], &
               [8,2], order=[2,1])
  y = [-1,-1,-1,-1,1,1,1,1]
  kernel = rbfdot(0.5_dp)

  call ksvm(x, y, kernel, model, cost=10.0_dp)
  call predict_kernel_model(model, x, scores, status, predicted)
  print '(a,8(1x,i0))', 'Predicted classes:', predicted

  call specc(x, 2, kernel, clusters)
  print '(a,8(1x,i0))', 'Spectral clusters:', clusters%labels
end program demo_kernlab
