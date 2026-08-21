program test_aalen
  use relsurv
  implicit none
  real(dp)::start(3),stop(3),x(3,0),times(1)
  integer::status(3)
  type(aalen_result)::fit
  start=0.0_dp; stop=[1.0_dp,2.0_dp,2.0_dp]; status=[1,0,0]; times=1.0_dp
  call aalen_fit(start,stop,status,x,times,fit,variance=.true.)
  if(abs(fit%increment(1,1)-1.0_dp/3.0_dp)>1.0e-12_dp) then
    print *, 'FAIL aalen ',fit%increment(1,1); error stop 1
  end if
  print *, 'test_aalen: PASS'
end program
