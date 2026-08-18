! SPDX-License-Identifier: GPL-2.0-or-later
program test_generators
  use runuran
  implicit none
  type(rng_state)::rng
  type(continuous_distribution)::d
  type(unuran_generator)::g
  real(dp),allocatable::x(:)
  real(dp)::m,s,c1,c2
  integer::fails
  fails=0;call rng_seed(rng,123456_i8);allocate(x(5000))
  d=udnorm();g=ars_new(d);call g%sample_n(rng,x)
  m=sum(x)/size(x);s=sqrt(sum((x-m)**2)/(size(x)-1))
  if(abs(m)>0.08_dp.or.abs(s-1.0_dp)>0.08_dp)then;print *,'ARS moments',m,s;fails=fails+1;end if
  d=udnorm(lb=-1.0_dp,ub=2.0_dp)
  g=pinv_new(d)
  if(.not.unuran_is_inversion(g))fails=fails+1
  if(abs(ud(g,0.0_dp)-d%pdf(0.0_dp))>1e-14_dp)fails=fails+1
  if(abs(up(g,0.0_dp)-d%cdf(0.0_dp))>1e-14_dp)fails=fails+1
  call g%sample_n(rng,x)
  if(any(x < -1.0_dp).or.any(x > 2.0_dp))then;print *,'truncation failure';fails=fails+1;end if
  c1=d%cdf(-1.0_dp)
  c2=d%cdf(2.0_dp)
  if(abs(c1)>1e-14_dp .or. abs(c2-1.0_dp)>1e-14_dp)then
    fails=fails+1
  end if
  if(fails==0)then;print *,'test_generators: PASS';else;error stop 1;end if
end program
