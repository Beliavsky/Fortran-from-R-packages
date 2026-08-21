program test_rsadd_residuals
  use relsurv, only : dp, rsadd_residual_result, rsadd_schoenfeld_residuals
  implicit none
  type(rsadd_residual_result) :: r
  real(dp) :: start(3), stop(3), z(3,1), beta(1), lp(3,1), lambda0(2), fup(2,1)
  integer :: status(3)
  start=0.0_dp
  stop=[1.0_dp,2.0_dp,3.0_dp]
  status=[1,1,0]
  z(:,1)=[0.0_dp,1.0_dp,2.0_dp]
  beta=0.0_dp;lp=0.0_dp;lambda0=1.0_dp;fup=0.0_dp
  call rsadd_schoenfeld_residuals(start,stop,status,z,beta,lp,lambda0,fup,r)
  call check(abs(r%residual(1,1)+1.0_dp)<1.0e-12_dp,'first residual')
  call check(abs(r%residual(2,1)+0.5_dp)<1.0e-12_dp,'second residual')
  call check(abs(r%var_event(1,1,1)-2.0_dp/3.0_dp)<1.0e-12_dp,'first variance')
  call check(abs(r%var_event(1,1,2)-0.25_dp)<1.0e-12_dp,'second variance')
  call check(abs(r%var_sum(1,1)-11.0_dp/12.0_dp)<1.0e-12_dp,'variance sum')
  call check(maxval(abs(r%kvar_event-r%var_event))<1.0e-12_dp,'partial variance no population')

  stop=[1.0_dp,1.0_dp,2.0_dp];status=[1,1,0];z(:,1)=[0.0_dp,2.0_dp,1.0_dp]
  call rsadd_schoenfeld_residuals(start,stop,status,z,beta,lp,lambda0,fup,r)
  call check(maxval(abs(r%residual(:,1)-[-1.0_dp,1.0_dp]))<1.0e-12_dp,'tied residuals')
  call check(maxval(abs(r%var_event(1,1,:)-2.0_dp/3.0_dp))<1.0e-12_dp,'tied variance averaging')
  print '(a)','test_rsadd_residuals: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then
      print '(a,1x,a)','FAIL',msg
      error stop 1
    end if
  end subroutine check
end program test_rsadd_residuals
