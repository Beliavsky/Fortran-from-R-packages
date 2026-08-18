! SPDX-License-Identifier: GPL-3.0-only
program test_generalized_reduction
  use poisson_binomial, only : dp, gpb_table, dpbinom, dgpbinom
  implicit none
  real(dp),parameter::p(5)=[0.11_dp,0.23_dp,0.48_dp,0.72_dp,0.91_dp]
  integer,parameter::vp(5)=[1,1,1,1,1],vq(5)=[0,0,0,0,0]
  real(dp),allocatable::a(:)
  type(gpb_table)::b
  a=dpbinom(p,"Characteristic")
  b=dgpbinom(p,vp,vq,"Characteristic")
  call chk(b%lower==0 .and. b%upper==5,"support")
  call chk(maxval(abs(a-b%values))<3e-13_dp,"ordinary reduction")
  print '(a)', 'test_generalized_reduction: PASS'
contains
  subroutine chk(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)',trim(msg); error stop 1; end if
  end subroutine chk
end program test_generalized_reduction
