! SPDX-License-Identifier: GPL-3.0-only
program test_tails_logs
  use poisson_binomial, only : dp, ppbinom_at, qpbinom, pgpbinom_at, qgpbinom
  implicit none
  real(dp),parameter::p(4)=[0.1_dp,0.3_dp,0.6_dp,0.8_dp]
  real(dp),parameter::gp(3)=[0.2_dp,0.7_dp,0.4_dp]
  integer,parameter::vp(3)=[3,-1,5],vq(3)=[0,2,1]
  real(dp)::v
  v=ppbinom_at(1,p,"Convolve",lower_tail=.false.)
  call chk(abs(v-0.6452_dp)<2e-14_dp,"ordinary upper tail")
  v=ppbinom_at(1,p,"Convolve",log_p=.true.)
  call chk(abs(v-log(0.3548_dp))<2e-14_dp,"ordinary log cdf")
  call chk(qpbinom(0.2_dp,p,"Convolve",lower_tail=.false.)==2,"ordinary upper q")
  v=pgpbinom_at(4,gp,vp,vq,"Convolve",lower_tail=.false.)
  call chk(abs(v-0.212_dp)<2e-14_dp,"generalized upper tail")
  v=pgpbinom_at(4,gp,vp,vq,"Convolve",log_p=.true.)
  call chk(abs(v-log(0.788_dp))<2e-14_dp,"generalized log cdf")
  call chk(qgpbinom(0.2_dp,gp,vp,vq,"Convolve",lower_tail=.false.)==6, &
           "generalized upper q")
  print '(a)', 'test_tails_logs: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_tails_logs
