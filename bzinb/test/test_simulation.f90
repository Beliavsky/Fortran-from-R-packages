program test_simulation
  use bzinb
  use test_support
  implicit none
  integer,parameter::n=50000
  integer::x(n),y(n),fail
  real(dp)::mx,my,p(4)
  fail=0;call set_bzinb_seed(12345)
  call rbp_sample(n,0.7_dp,1.2_dp,0.8_dp,x,y)
  mx=sum(real(x,dp))/n;my=sum(real(y,dp))/n
  call assert_close('rbp mean x',mx,1.9_dp,0.035_dp,fail)
  call assert_close('rbp mean y',my,1.5_dp,0.035_dp,fail)
  p=[0.6_dp,0.15_dp,0.15_dp,0.10_dp]
  call rbzinb_sample(n,1.5_dp,0.8_dp,1.1_dp,0.7_dp,1.0_dp,p,x,y)
  mx=sum(real(x,dp))/n;my=sum(real(y,dp))/n
  call assert_true('rbzinb finite means',mx>0.1_dp.and.my>0.1_dp,fail)
  if(fail==0)then;print *,'test_simulation: PASS';else;error stop 1;end if
end program
