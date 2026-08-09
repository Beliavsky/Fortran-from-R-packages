program test_ed
  use smoof_kinds, only : dp
  use smoof_ed
  implicit none
  real(dp)::f(3),x(5),theta(2)
  x=[0.1_dp,0.2_dp,0.0_dp,0.0_dp,0.0_dp]
  theta=[0.3_dp,0.4_dp]
  call ed1(x,3,2.0_dp,theta,f)
  if(any(abs(f)>=huge(1.0_dp)))error stop 'ed1 nonfinite'
  call ed2(x,3,2.0_dp,theta,f)
  if(any(abs(f)>=huge(1.0_dp)))error stop 'ed2 nonfinite'
  print *,'test_ed: PASS'
end program test_ed
