program test_multi
  use smoof_kinds, only : dp
  use smoof_multi
  implicit none
  real(dp)::f2(2),f3(3),x7(7),z6(6)
  integer::fails
  fails=0
  x7=0.5_dp
  call dtlz2(x7,3,f3)
  call vcheck('dtlz2',f3,[0.5_dp,0.5_dp,sqrt(0.5_dp)],1.0e-12_dp,fails)
  call zdt1([0.25_dp,0.0_dp,0.0_dp],f2)
  call vcheck('zdt1',f2,[0.25_dp,0.5_dp],1.0e-12_dp,fails)
  call mop1([1.0_dp,0.0_dp],f2)
  call vcheck('mop1',f2,[1.0_dp,1.0_dp],1.0e-14_dp,fails)
  call bk1([1.0_dp,2.0_dp],f2)
  call vcheck('bk1',f2,[3.0_dp,25.0_dp],1.0e-14_dp,fails)
  z6=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
  call wfg1(z6,3,2,f3)
  call vcheck('wfg1',f3,[2.886792851925874_dp,0.973268463057909_dp, &
    0.974904813720709_dp],2.0e-12_dp,fails)
  call wfg9(z6,3,2,f3)
  if (any(abs(f3)>=huge(1.0_dp))) fails=fails+1
  if(fails/=0) error stop 'test_multi failed'
  print *, 'test_multi: PASS'
contains
  subroutine vcheck(name,a,b,tol,fails)
    character(*),intent(in)::name
    real(dp),intent(in)::a(:),b(:),tol
    integer,intent(inout)::fails
    if(size(a)/=size(b).or.any(abs(a-b)>tol))then
      print *,'FAIL ',trim(name),a,b;fails=fails+1
    end if
  end subroutine vcheck
end program test_multi
