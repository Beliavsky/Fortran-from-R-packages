program test_ajfit
  use flexsurv_kinds, only : dp
  use flexsurv_ajfit, only : ajfit_result, ajfit_transition_data
  implicit none
  type(ajfit_result)::aj
  real(dp)::time(6)
  integer::status(6),from(6),to(6),fails
  fails=0
  time=[1.0_dp,2.0_dp,3.0_dp,1.0_dp,2.0_dp,3.0_dp]
  status=[1,0,0,0,1,0]
  from=1;to=[2,2,2,3,3,3]
  aj=ajfit_transition_data(time,status,from,to,3,1)
  if(size(aj%time)/=3)then;print *,'FAIL size';fails=fails+1;end if
  call chk(aj%probability(1,2),2.0_dp/3.0_dp,'p1 t1',fails)
  call chk(aj%probability(2,2),1.0_dp/3.0_dp,'p2 t1',fails)
  call chk(aj%probability(1,3),1.0_dp/3.0_dp,'p1 t2',fails)
  call chk(aj%probability(2,3),1.0_dp/3.0_dp,'p2 t2',fails)
  call chk(aj%probability(3,3),1.0_dp/3.0_dp,'p3 t2',fails)
  if(fails==0)then;print *,'test_ajfit: PASS';else;print *,'test_ajfit: FAIL',fails;error stop 1;end if
contains
  subroutine chk(a,b,name,fails)
    real(dp),intent(in)::a,b
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    if(abs(a-b)>1.0e-12_dp)then;print *,'FAIL ',trim(name),a,b;fails=fails+1;end if
  end subroutine chk
end program test_ajfit
