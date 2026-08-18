! SPDX-License-Identifier: GPL-3.0-only
program test_generalized_gcd_weights
  use poisson_binomial, only : dp, gpb_table, dgpbinom, dgpbinom_at, qgpbinom
  implicit none
  real(dp), parameter :: p(2) = [0.25_dp,0.6_dp]
  integer, parameter :: vp(2) = [2,4], vq(2) = [0,0], w(2) = [2,1]
  type(gpb_table) :: a, b, c
  a=dgpbinom(p,vp,vq,"Convolve",w)
  b=dgpbinom(p,vp,vq,"DivideFFT",w)
  c=dgpbinom(p,vp,vq,"Characteristic",w)
  call chk(a%lower==0 .and. a%upper==8,"support")
  call chk(abs(dgpbinom_at(0,p,vp,vq,"Convolve",w)-0.225_dp)<2e-14_dp,"p0")
  call chk(abs(dgpbinom_at(2,p,vp,vq,"Convolve",w)-0.15_dp)<2e-14_dp,"p2")
  call chk(abs(dgpbinom_at(4,p,vp,vq,"Convolve",w)-0.3625_dp)<2e-14_dp,"p4")
  call chk(abs(dgpbinom_at(6,p,vp,vq,"Convolve",w)-0.225_dp)<2e-14_dp,"p6")
  call chk(abs(dgpbinom_at(8,p,vp,vq,"Convolve",w)-0.0375_dp)<2e-14_dp,"p8")
  call chk(abs(dgpbinom_at(1,p,vp,vq,"Convolve",w))<2e-14_dp,"gcd gap")
  call chk(maxval(abs(a%values-b%values))<3e-13_dp,"divide agreement")
  call chk(maxval(abs(a%values-c%values))<3e-13_dp,"characteristic agreement")
  call chk(qgpbinom(0.70_dp,p,vp,vq,"Convolve",w)==4,"weighted quantile")
  print '(a)', 'test_generalized_gcd_weights: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_generalized_gcd_weights
