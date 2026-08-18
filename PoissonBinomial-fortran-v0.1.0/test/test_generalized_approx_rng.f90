! SPDX-License-Identifier: GPL-3.0-only
program test_generalized_approx_rng
  use poisson_binomial, only : dp, gpb_table, pgpbinom, qgpbinom, rgpbinom
  implicit none
  real(dp),parameter::p(3)=[0.2_dp,0.7_dp,0.4_dp]
  integer,parameter::vp(3)=[3,-1,5],vq(3)=[0,2,1]
  type(gpb_table)::a
  integer,allocatable::x(:)
  real(dp)::avg
  a=pgpbinom(p,vp,vq,"Normal",lower_tail=.true.)
  call chk(abs(a%values(6)-0.8149527658498652_dp)<3e-13_dp,"normal")
  a=pgpbinom(p,vp,vq,"RefinedNormal",lower_tail=.true.)
  call chk(abs(a%values(6)-0.8185678474026189_dp)<3e-13_dp,"refined")
  call chk(qgpbinom(0.75_dp,p,vp,vq,"Convolve")==4,"quantile")
  call rgpbinom(20000,p,vp,vq,x,"Convolve",generator="Bernoulli")
  avg=sum(real(x,dp))/real(size(x),dp)
  call chk(abs(avg-3.1_dp)<0.06_dp,"rng mean")
  call chk(minval(x)>=0 .and. maxval(x)<=10,"rng support")
  print '(a)', 'test_generalized_approx_rng: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_generalized_approx_rng
