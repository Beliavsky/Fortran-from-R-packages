program test_km
  use survival
  implicit none
  type(survfit_result) :: fit
  real(dp) :: time(4)=[1.0_dp,2.0_dp,2.0_dp,3.0_dp]
  integer :: status(4)=[1,1,0,1]
  call kaplan_meier(time,status,fit)
  if(size(fit%time)/=3) error stop 'km times'
  if(abs(fit%survival(1)-0.75_dp)>1e-12_dp) error stop 'km s1'
  if(abs(fit%survival(2)-0.5_dp)>1e-12_dp) error stop 'km s2'
  if(abs(fit%cumhaz(2)-(0.25_dp+1.0_dp/3.0_dp))>1e-12_dp) error stop 'km h2'
  print *, 'test_km PASS'
end program
