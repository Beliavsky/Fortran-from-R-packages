! SPDX-License-Identifier: GPL-3.0-only
program test_generalized_exact
  use poisson_binomial, only : dp, gpb_table, dgpbinom, dgpbinom_at, dgpbinom_values, &
                               pgpbinom_at, pgpbinom_values
  implicit none
  real(dp),parameter::p(3)=[0.2_dp,0.7_dp,0.4_dp]
  integer,parameter::vp(3)=[3,-1,5],vq(3)=[0,2,1]
  type(gpb_table)::a,b,c
  real(dp),allocatable::v(:)
  a=dgpbinom(p,vp,vq,"Convolve")
  b=dgpbinom(p,vp,vq,"DivideFFT")
  c=dgpbinom(p,vp,vq,"Characteristic")
  call chk(a%lower==0 .and. a%upper==10,"support")
  call chk(abs(dgpbinom_at(0,p,vp,vq,"Convolve")-0.336_dp)<2e-14_dp,"p0")
  call chk(abs(dgpbinom_at(3,p,vp,vq,"Convolve")-0.228_dp)<2e-14_dp,"p3")
  call chk(abs(dgpbinom_at(4,p,vp,vq,"Convolve")-0.224_dp)<2e-14_dp,"p4")
  call chk(abs(dgpbinom_at(6,p,vp,vq,"Convolve")-0.036_dp)<2e-14_dp,"p6")
  call chk(abs(dgpbinom_at(7,p,vp,vq,"Convolve")-0.152_dp)<2e-14_dp,"p7")
  call chk(abs(dgpbinom_at(10,p,vp,vq,"Convolve")-0.024_dp)<2e-14_dp,"p10")
  call chk(maxval(abs(a%values-b%values))<2e-13_dp,"divide agreement")
  call chk(maxval(abs(a%values-c%values))<3e-13_dp,"characteristic agreement")
  call chk(abs(pgpbinom_at(4,p,vp,vq,"Convolve")-0.788_dp)<2e-14_dp,"cdf")
  v=dgpbinom_values([0,3,4,7],p,vp,vq,"Convolve")
  call chk(maxval(abs(v-[0.336_dp,0.228_dp,0.224_dp,0.152_dp]))<2e-14_dp,"vector density")
  v=pgpbinom_values([0,3,4],p,vp,vq,"Convolve")
  call chk(maxval(abs(v-[0.336_dp,0.564_dp,0.788_dp]))<2e-14_dp,"vector cdf")
  print '(a)', 'test_generalized_exact: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_generalized_exact
