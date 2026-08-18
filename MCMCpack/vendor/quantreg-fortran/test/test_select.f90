program test_select
  use quantreg, only : dp, kuantiles, qselect
  implicit none
  real(dp) :: x(5), p(3), q(3)
  x = [5.0_dp,1.0_dp,4.0_dp,2.0_dp,3.0_dp]
  p = [0.25_dp,0.5_dp,0.75_dp]
  call kuantiles(x,p,q,7)
  if (maxval(abs(q-[2.0_dp,3.0_dp,4.0_dp])) > 1.0e-12_dp) error stop 'type 7'
  if (abs(qselect(x,0.5_dp)-3.0_dp) > 1.0e-12_dp) error stop 'qselect'
  print *, 'test_select: PASS'
end program
