program test_cec2019
  use smoof_kinds, only : dp
  use smoof_cec2019
  implicit none
  real(dp)::f2(2),f3(3),x3(3)
  integer::fails
  fails=0
  call mmf1([2.0_dp,0.0_dp],f2)
  call check('mmf1',f2,[0.0_dp,1.0_dp],1.0e-13_dp,fails)
  call omni_test([0.0_dp,0.0_dp],f2)
  call check('omni',f2,[0.0_dp,2.0_dp],1.0e-13_dp,fails)
  x3=[0.2_dp,0.3_dp,0.4_dp]
  call mmf14(x3,3,2,f3)
  if(any(abs(f3)>=huge(1.0_dp))) fails=fails+1
  if(fails/=0) error stop 'test_cec2019 failed'
  print *,'test_cec2019: PASS'
contains
  subroutine check(name,a,b,tol,fails)
    character(*),intent(in)::name;real(dp),intent(in)::a(:),b(:),tol
    integer,intent(inout)::fails
    if(any(abs(a-b)>tol))then;print *,'FAIL ',trim(name),a,b;fails=fails+1;end if
  end subroutine check
end program test_cec2019
