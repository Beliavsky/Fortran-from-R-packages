! SPDX-License-Identifier: GPL-3.0-only
program test_random
  use, intrinsic :: iso_fortran_env, only: int64
  use hdshop, only: dp, random_covariance_matrix, symmetric_eigen
  implicit none
  real(dp),allocatable::a(:,:),b(:,:),ev(:),v(:,:)
  logical::ok
  a=random_covariance_matrix(5,[0.2_dp,0.4_dp,0.8_dp,1.6_dp,3.2_dp],1234_int64)
  b=random_covariance_matrix(5,[0.2_dp,0.4_dp,0.8_dp,1.6_dp,3.2_dp],1234_int64)
  if(maxval(abs(a-b))>0.0_dp)error stop 1
  call symmetric_eigen(a,ev,v,ok)
  if(.not.ok .or. maxval(abs(ev-[0.2_dp,0.4_dp,0.8_dp,1.6_dp,3.2_dp]))>1.0e-10_dp)error stop 1
  print '(a)', 'test_random: PASS'
end program test_random
