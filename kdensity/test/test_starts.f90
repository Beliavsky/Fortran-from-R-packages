program test_starts
  use kdensity
  implicit none
  type(kd_start)::st
  real(dp),allocatable::p(:)
  real(dp)::x(5)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  integer::s
  st=get_start('normal',s);call st%estimator(x,p,s);call check(s==0.and.size(p)==2,'normal fit');call check(abs(p(1)-3.0_dp)<1e-14_dp,'normal mean')
  st=get_start('exponential',s);call st%estimator(x,p,s);call check(abs(p(1)-1.0_dp/3.0_dp)<1e-14_dp,'exp rate')
  st=get_start('gamma',s);call st%estimator(x,p,s);call check(all(p>0),'gamma fit')
  st=get_start('weibull',s);call st%estimator(x,p,s);call check(all(p>0),'weibull fit')
  st=get_start('pareto',s);call st%estimator(x,p,s);call check(p(1)>0,'pareto fit')
  print *, 'test_starts: PASS'
contains
  subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg;if(.not.ok)error stop msg;end subroutine
end program
