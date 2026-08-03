! SPDX-License-Identifier: GPL-2.0-only
program svm_gp
  use kernlab
  implicit none
  real(dp) :: x(8,2), yr(8)
  integer :: yc(8), status
  type(kernel_spec) :: kernel
  type(kernel_model) :: svm_model, gp_model
  real(dp), allocatable :: scores(:,:), variance(:)
  integer, allocatable :: classes(:)

  x = reshape([ -2.0_dp,-2.0_dp, -1.8_dp,-2.2_dp, -2.2_dp,-1.8_dp, -1.9_dp,-1.7_dp, &
                 2.0_dp, 2.0_dp,  1.8_dp, 2.2_dp,  2.2_dp, 1.8_dp,  1.9_dp, 1.7_dp ], &
               [8,2], order=[2,1])
  yc = [-1,-1,-1,-1,1,1,1,1]
  yr = 2.0_dp*x(:,1)-x(:,2)
  kernel = rbfdot(0.5_dp)

  call ksvm(x, yc, kernel, svm_model, cost=10.0_dp)
  call predict_kernel_model(svm_model, x, scores, status, classes)
  print '(a,8(1x,i0))', 'SVM classes:', classes

  call gausspr(x, yr, kernel, gp_model, var=1.0e-4_dp)
  call gausspr_predict_variance(gp_model, x, scores, variance, status)
  print '(a)', 'Gaussian-process means and variances:'
  print '(2f12.6)', [scores(:,1),variance]
end program svm_gp
