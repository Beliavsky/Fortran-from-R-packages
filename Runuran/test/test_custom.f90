! SPDX-License-Identifier: GPL-2.0-or-later
program test_custom
  use runuran
  implicit none
  type(continuous_distribution)::d
  type(unuran_generator)::g
  type(rng_state)::rng
  real(dp)::x(2000),m
  d=ud_continuous(epdf,ecdf,0.0_dp,20.0_dp,name='custom exponential')
  if(abs(d%cdf(log(2.0_dp))-0.5_dp)>2e-12_dp)then
    print *,'custom cdf failed',d%cdf(log(2.0_dp));error stop 1
  end if
  g=pinv_new(d);call rng_seed(rng,31415_i8);call g%sample_n(rng,x);m=sum(x)/size(x)
  if(abs(m-1.0_dp)>0.08_dp)then;print *,'custom sample mean',m;error stop 1;end if
  print *,'test_custom: PASS'
contains
  real(dp) function epdf(x,p) result(y)
    real(dp),intent(in)::x,p(:)
    if(size(p)<0) stop
    if(x<0.0_dp)then;y=0.0_dp;else;y=exp(-x);end if
  end function
  real(dp) function ecdf(x,p) result(y)
    real(dp),intent(in)::x,p(:)
    if(size(p)<0) stop
    if(x<=0.0_dp)then;y=0.0_dp;else;y=1.0_dp-exp(-x);end if
  end function
end program
