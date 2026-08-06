program test_fit_bounded
  use kdensity
  implicit none
  real(dp)::xb(8)=[0.08_dp,0.12_dp,0.18_dp,0.25_dp,0.32_dp,0.45_dp,0.60_dp,0.72_dp]
  real(dp)::xg(7)=[0.2_dp,0.4_dp,0.7_dp,1.0_dp,1.4_dp,2.0_dp,2.8_dp]
  type(kdensity_fit)::fit
  type(kdensity_options)::opt
  real(dp)::area
  integer::s
  opt%kernel='beta';opt%start='uniform';opt%bandwidth='HS';opt%support=[0.0_dp,1.0_dp];opt%support_supplied=.true.
  fit=fit_kdensity(xb,opt);call check(fit%status==0,trim(fit%message));area=adaptive_integral(fb,0.0_dp,1.0_dp,1e-6_dp,s);call check(abs(area-1)<2e-4_dp,'beta normalized')
  opt%kernel='gamma';opt%start='gamma';opt%bandwidth='nrd0';opt%support=[0.0_dp,huge(1.0_dp)];opt%support_supplied=.true.
  fit=fit_kdensity(xg,opt);call check(fit%status==0,trim(fit%message));call check(fit%pdf(1.0_dp)>0,'gamma density')
  print *, 'test_fit_bounded: PASS'
contains
  function fb(y) result(v);real(dp),intent(in)::y;real(dp)::v;v=fit%pdf(y);end function
  subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg;if(.not.ok)error stop msg;end subroutine
end program
