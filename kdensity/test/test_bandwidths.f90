program test_bandwidths
  use kdensity
  implicit none
  real(dp)::x(10)=[0.10_dp,0.15_dp,0.20_dp,0.25_dp,0.30_dp,0.35_dp,0.40_dp,0.50_dp,0.65_dp,0.80_dp]
  real(dp)::h
  integer::s
  logical::fallback
  call check(bandwidth_nrd0(x)>0,'nrd0')
  call check(bandwidth_nrd(x)>0,'nrd')
  call check(bandwidth_jh(x)>0,'jh')
  call check(bandwidth_rhe(x)>0,'rhe')
  h=bandwidth_hs(x,s,fallback);call check(s==0.and.h>0,'hs')
  print *, 'test_bandwidths: PASS'
contains
  subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg;if(.not.ok)error stop msg;end subroutine
end program
