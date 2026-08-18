! SPDX-License-Identifier: GPL-3.0-only
program test_weights_quantile_rng
  use poisson_binomial, only : dp, dpbinom, dpbinom_values, ppbinom_at, ppbinom_values, &
                               qpbinom, qpbinom_values, rpbinom
  implicit none
  real(dp),parameter::p(3)=[0.0_dp,0.25_dp,1.0_dp]
  integer,parameter::w(3)=[2,2,1]
  real(dp),allocatable::pmf(:),v(:)
  integer,allocatable::x(:),qv(:)
  real(dp)::avg
  pmf=dpbinom(p,"Convolve",w)
  call chk(size(pmf)==6,"weighted support")
  call chk(abs(pmf(2)-0.5625_dp)<1e-14_dp,"weighted p1")
  call chk(abs(pmf(3)-0.375_dp)<1e-14_dp,"weighted p2")
  call chk(abs(pmf(4)-0.0625_dp)<1e-14_dp,"weighted p3")
  call chk(abs(ppbinom_at(2,p,"Convolve",w)-0.9375_dp)<1e-14_dp,"cdf")
  call chk(qpbinom(0.90_dp,p,"Convolve",w)==2,"quantile")
  v=dpbinom_values([1,2,3],p,"Convolve",w)
  call chk(maxval(abs(v-[0.5625_dp,0.375_dp,0.0625_dp]))<1e-14_dp,"vector density")
  v=ppbinom_values([0,1,2],p,"Convolve",w)
  call chk(maxval(abs(v-[0.0_dp,0.5625_dp,0.9375_dp]))<1e-14_dp,"vector cdf")
  qv=qpbinom_values([0.1_dp,0.9_dp],p,"Convolve",w)
  call chk(all(qv==[1,2]),"vector quantile")
  call rpbinom(20000,p,x,"Convolve",w,"Bernoulli")
  avg=sum(real(x,dp))/real(size(x),dp)
  call chk(abs(avg-1.5_dp)<0.035_dp,"rng mean")
  call chk(minval(x)>=1 .and. maxval(x)<=3,"rng support")
  print '(a)', 'test_weights_quantile_rng: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_weights_quantile_rng
