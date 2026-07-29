! SPDX-License-Identifier: MIT
program test_utils
  use bayesianou
  implicit none
  real(dp) :: x(5,2),factor(5),load(2),beta(2),res(5)
  type(zscore_result) :: z
  type(ou_level_spec_type) :: spec
  integer :: status
  x(:,1)=[1,2,3,4,5];x(:,2)=[2,4,6,8,10]
  call zscore_train(x,4,z)
  call assert_close(sum(z%mz(1:4,1)),0.0_dp,1e-12_dp)
  call assert_close(z%sd(1),sqrt(5.0_dp/3.0_dp),1e-12_dp)
  call compute_common_factor(z%mz,4,.true.,factor,load,status)
  if(status/=status_ok)error stop 'PCA failed'
  if(abs(abs(load(1))-1.0_dp/sqrt(2.0_dp))>1e-6_dp) error stop 'bad PCA loading'
  call orthogonalize_series(factor,factor,4,res,beta,status)
  if(maxval(abs(res))>1e-7_dp)error stop 'orthogonalization failed'
  spec=ou_level_spec(level_canonical)
  if(.not.spec%level1%cubic.or.spec%level2%cubic)error stop 'bad level spec'
  print *, 'test_utils: PASS'
contains
  subroutine assert_close(a,b,tol)
    real(dp),intent(in)::a,b,tol
    if(abs(a-b)>tol)error stop 'assert_close failed'
  end subroutine assert_close
end program test_utils
