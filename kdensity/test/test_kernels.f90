program test_kernels
  use kdensity
  implicit none
  type(kd_kernel) :: k
  integer :: s
  real(dp) :: v
  k=get_kernel('gaussian',s); call check(s==0,'gaussian lookup')
  call check(abs(k%evaluate(0.0_dp,0.0_dp,1.0_dp)-normal_pdf(0.0_dp))<1e-14_dp,'gaussian value')
  k=get_kernel('epanechnikov',s);call check(abs(k%evaluate(0.0_dp,0.0_dp,1.0_dp)-0.75_dp)<1e-14_dp,'epanechnikov')
  k=get_kernel('gamma',s);v=k%evaluate(0.5_dp,0.5_dp,0.1_dp);call check(v>0,'gamma kernel')
  k=get_kernel('beta',s);v=k%evaluate(0.2_dp,0.3_dp,0.05_dp);call check(v>0,'beta kernel')
  k=get_kernel('gcopula',s);v=k%evaluate(0.3_dp,0.4_dp,0.2_dp);call check(v>0,'copula kernel')
  print *, 'test_kernels: PASS'
contains
  subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg;if(.not.ok)error stop msg;end subroutine
end program
