program test_cec09
  use smoof_kinds, only : dp
  use smoof_cec09
  implicit none
  real(dp)::x(5),f2(2),f3(3)
  integer::fails
  fails=0
  x=[0.25_dp,0.1_dp,-0.2_dp,0.3_dp,-0.1_dp]
  call uf1(x,f2)
  call check('uf1',f2,[1.719098300562505_dp,0.926393202250021_dp],2.0e-12_dp,fails)
  x=[0.2_dp,0.4_dp,-0.2_dp,0.1_dp,-0.3_dp]
  call uf8(x,f3)
  if(any(abs(f3)>=huge(1.0_dp)))fails=fails+1
  call uf10(x,f3)
  if(any(abs(f3)>=huge(1.0_dp)))fails=fails+1
  if(fails/=0)error stop 'test_cec09 failed'
  print *,'test_cec09: PASS'
contains
  subroutine check(name,a,b,tol,fails)
    character(*),intent(in)::name;real(dp),intent(in)::a(:),b(:),tol
    integer,intent(inout)::fails
    if(any(abs(a-b)>tol))then;print *,'FAIL ',trim(name),a,b;fails=fails+1;end if
  end subroutine check
end program test_cec09
