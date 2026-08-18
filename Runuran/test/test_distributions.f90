! SPDX-License-Identifier: GPL-2.0-or-later
program test_distributions
  use runuran
  implicit none
  type(continuous_distribution)::d
  integer::fails
  fails=0
  d=udnorm()
  call chk(d%pdf(0.0_dp),0.3989422804014327_dp,2e-14_dp,'normal pdf')
  call chk(d%cdf(0.0_dp),0.5_dp,2e-15_dp,'normal cdf')
  call chk(d%quantile(0.975_dp),1.959963984540054_dp,2e-9_dp,'normal quantile')
  d=udgamma(2.0_dp,3.0_dp)
  call chk(d%cdf(6.0_dp),1.0_dp-3.0_dp*exp(-2.0_dp),2e-13_dp,'gamma cdf')
  d=udbeta(2.0_dp,3.0_dp)
  call chk(d%cdf(0.5_dp),0.6875_dp,2e-13_dp,'beta cdf')
  d=udcauchy(2.0_dp,3.0_dp)
  call chk(d%quantile(0.5_dp),2.0_dp,2e-14_dp,'cauchy median')
  d=udpowerexp(2.0_dp)
  call chk(d%pdf(0.0_dp),1.0_dp/sqrt(pi),2e-13_dp,'power exponential')
  d=udburr(2.0_dp,3.0_dp)
  call chk(d%cdf(1.0_dp),0.75_dp,2e-13_dp,'Burr cdf')
  if(fails==0)then;print *,'test_distributions: PASS';else;error stop 1;end if
contains
  subroutine chk(x,e,t,n)
    real(dp),intent(in)::x,e,t;character(len=*),intent(in)::n
    if(abs(x-e)>t*max(1.0_dp,abs(e)))then
      print *,trim(n),' got=',x,' expected=',e;fails=fails+1
    end if
  end subroutine
end program
