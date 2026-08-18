! SPDX-License-Identifier: GPL-2.0-or-later
program test_mixture
  use runuran
  implicit none
  type(rng_state)::rng
  type(continuous_distribution)::c(2)
  type(unuran_generator)::g
  real(dp),allocatable::x(:)
  real(dp)::m
  call rng_seed(rng,777_i8);c(1)=udnorm(-2.0_dp,1.0_dp);c(2)=udnorm(3.0_dp,1.0_dp)
  g=mixt_new(c,[0.4_dp,0.6_dp]);allocate(x(10000));call g%sample_n(rng,x);m=sum(x)/size(x)
  if(abs(m-1.0_dp)>0.12_dp)then;print *,'mixture mean',m;error stop 1;end if
  print *,'test_mixture: PASS'
end program
