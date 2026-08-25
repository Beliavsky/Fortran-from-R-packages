program test_factorial
  use compound_cox
  implicit none
  real(dp) :: t(9)
  integer :: d(9), g(9)
  type(factorial_result) :: fr
  t=[1.0_dp,2.0_dp,4.0_dp,1.5_dp,3.0_dp,5.0_dp,2.2_dp,4.2_dp,6.0_dp]
  d=[1,1,0,1,0,1,1,0,1]
  g=[1,1,1,2,2,2,3,3,3]
  call surv_factorial(t,d,g,2.0_dp,fr,copula='clayton',nsim=100)
  if(size(fr%estimate)/=3) error stop 'factorial size'
  if(any(fr%estimate<0.0_dp) .or. any(fr%estimate>1.0_dp)) error stop 'factorial estimate range'
  if(fr%p_simu<0.0_dp .or. fr%p_simu>1.0_dp) error stop 'factorial simu p'
  if(fr%p_anal<0.0_dp .or. fr%p_anal>1.0_dp) error stop 'factorial anal p'
  print '(a)', 'test_factorial: PASS'
end program test_factorial
