! SPDX-License-Identifier: GPL-3.0-only
program test_sample_generators
  use poisson_binomial, only : dp, rpbinom, rgpbinom
  implicit none
  real(dp), parameter :: p(4)=[0.1_dp,0.3_dp,0.6_dp,0.8_dp]
  integer, parameter :: vp(4)=[2,4,-1,3],vq(4)=[0,1,2,3]
  integer, allocatable :: x(:)
  real(dp) :: avg, expected
  call rpbinom(15000,p,x,"Convolve",generator="Sample")
  avg=sum(real(x,dp))/real(size(x),dp)
  call chk(abs(avg-sum(p))<0.04_dp,"ordinary sample mean")
  call chk(minval(x)>=0 .and. maxval(x)<=4,"ordinary sample support")
  expected=sum(p*real(vp,dp)+(1.0_dp-p)*real(vq,dp))
  call rgpbinom(15000,p,vp,vq,x,"Convolve",generator="Sample")
  avg=sum(real(x,dp))/real(size(x),dp)
  call chk(abs(avg-expected)<0.08_dp,"generalized sample mean")
  call chk(minval(x)>=sum(min(vp,vq)) .and. maxval(x)<=sum(max(vp,vq)), &
           "generalized sample support")
  print '(a)', 'test_sample_generators: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_sample_generators
