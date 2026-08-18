! SPDX-License-Identifier: GPL-2.0-or-later
program test_discrete
  use runuran
  implicit none
  type(discrete_distribution)::d
  integer::fails
  fails=0
  d=udbinom(10,0.2_dp)
  call chk(d%pmf(3),0.201326592_dp,2e-13_dp,'binomial pmf')
  call chk(d%cdf(2),0.6777995264_dp,2e-13_dp,'binomial cdf')
  if(d%quantile(0.7_dp)/=3)then;print *,'binomial quantile failed';fails=fails+1;end if
  d=udpois(2.0_dp)
  call chk(d%pmf(0),exp(-2.0_dp),2e-13_dp,'poisson pmf')
  d=udgeom(0.25_dp)
  call chk(d%pmf(2),0.140625_dp,2e-14_dp,'geometric pmf')
  d=udlogarithmic(0.5_dp)
  call chk(d%pmf(1),-0.5_dp/log(0.5_dp),2e-13_dp,'logarithmic pmf')
  if(fails==0)then;print *,'test_discrete: PASS';else;error stop 1;end if
contains
  subroutine chk(x,e,t,n)
    real(dp),intent(in)::x,e,t;character(len=*),intent(in)::n
    if(abs(x-e)>t*max(1.0_dp,abs(e)))then;print *,trim(n),x,e;fails=fails+1;end if
  end subroutine
end program
