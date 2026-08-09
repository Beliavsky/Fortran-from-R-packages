program test_nk
  use smoof_kinds, only : dp
  use smoof_nk
  implicit none
  integer :: x(2),k(2),ls(2),links(2),vs(2)
  real(dp)::vals(8),f
  k=[1,1];ls=[1,2];links=[2,1];vs=[1,5]
  vals=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,10.0_dp,20.0_dp,30.0_dp,40.0_dp]
  x=[1,0]
  f=nk_evaluate_raw(x,k,ls,links,vs,vals)
  if(abs(f-11.0_dp)>1.0e-14_dp)then
    print *,'FAIL nk',f;error stop
  end if
  print *,'test_nk: PASS'
end program test_nk
