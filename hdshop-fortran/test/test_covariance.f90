! SPDX-License-Identifier: GPL-3.0-only
program test_covariance
  use hdshop, only: dp, sigma_sample_estimator, cov_shrink_bgp14, &
    inv_cov_shrink_bgp16, nonlin_shrink_lw, matrix_shrink_result, &
    inverse_matrix, symmetric_eigen
  implicit none
  real(dp) :: x(3,8), target(3,3)
  real(dp), allocatable :: s(:,:), inv(:,:), nl(:,:), ev(:), vec(:,:)
  type(matrix_shrink_result) :: cr, ir
  logical :: ok

  x = reshape([ &
    0.01_dp,0.005_dp,-0.002_dp, 0.02_dp,-0.01_dp,0.008_dp, &
   -0.01_dp,0.015_dp,0.012_dp, 0.03_dp,0.02_dp,-0.004_dp, &
    0.00_dp,-0.005_dp,0.016_dp, 0.015_dp,0.01_dp,0.006_dp, &
   -0.005_dp,0.012_dp,0.011_dp, 0.025_dp,0.018_dp,0.003_dp], [3,8])
  s=sigma_sample_estimator(x)
  call assert_close(s(1,1),2.102678571428571e-4_dp,1.0e-16_dp)
  call assert_close(s(2,3),-3.289285714285714e-5_dp,1.0e-16_dp)
  target=0.0_dp;target(1,1)=1.0_dp/3.0_dp;target(2,2)=target(1,1);target(3,3)=target(1,1)
  cr=cov_shrink_bgp14(8,target,s)
  if(.not.cr%ok)error stop 1
  call assert_close(cr%alpha,0.358943665705116_dp,2.0e-13_dp)
  call assert_close(cr%beta,2.40327440753336e-4_dp,2.0e-16_dp)
  call inverse_matrix(s,inv,ok);if(.not.ok)error stop 1
  ir=inv_cov_shrink_bgp16(8,target,inv)
  if(.not.ir%ok)error stop 1
  call assert_close(ir%alpha,0.313024386269686_dp,2.0e-12_dp)
  nl=nonlin_shrink_lw(x)
  call symmetric_eigen(nl,ev,vec,ok)
  if(.not.ok .or. minval(ev)<=0.0_dp)error stop 1
  if(maxval(abs(nl-transpose(nl)))>1.0e-13_dp)error stop 1
  print '(a)', 'test_covariance: PASS'
contains
  subroutine assert_close(a,b,tol)
    real(dp),intent(in)::a,b,tol
    if(abs(a-b)>tol)then
      print *, 'mismatch:',a,b,abs(a-b);error stop 1
    end if
  end subroutine assert_close
end program test_covariance
