program test_fit_symmetric
  use kdensity
  implicit none
  real(dp)::x(7)=[-1.2_dp,-0.7_dp,-0.2_dp,0.0_dp,0.3_dp,0.8_dp,1.1_dp]
  type(kdensity_fit)::fit
  type(kdensity_options)::opt
  real(dp)::area
  integer::s
  opt%kernel='gaussian';opt%start='normal';opt%bandwidth='nrd0';opt%normalized=.true.
  fit=fit_kdensity(x,opt);call check(fit%status==0,trim(fit%message));call check(fit%pdf(0.0_dp)>0,'positive density')
  area=adaptive_integral(eval,-8.0_dp,8.0_dp,1e-7_dp,s);call check(abs(area-1)<2e-4_dp,'normalization')
  opt%bw=huge(1.0_dp);fit=fit_kdensity(x,opt);call check(abs(fit%pdf(0.2_dp)-fit%start_pdf(0.2_dp))<1e-14_dp,'infinite bandwidth')
  print *, 'test_fit_symmetric: PASS'
contains
  function eval(y) result(v);real(dp),intent(in)::y;real(dp)::v;v=fit%pdf(y);end function
  subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg;if(.not.ok)error stop msg;end subroutine
end program
