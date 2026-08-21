program test_dde_facade
  use desolve, only : dp,ode_result,history_buffer,dede_rk4,ode
  implicit none
  type(history_buffer)::hist
  type(ode_result)::d,s
  real(dp)::tt(3)
  call hist%init(1,1000);call hist%append(-1.0_dp,[1.0_dp],[0.0_dp]);call hist%append(0.0_dp,[1.0_dp],[0.0_dp])
  tt=[0.0_dp,0.25_dp,0.5_dp]
  d=dede_rk4(dr,[1.0_dp],tt,h=0.005_dp,initial_history=hist)
  if(.not.d%ok())error stop 'dde status'
  if(abs(d%y(1,3)-0.5_dp)>2e-4_dp)then;print *,d%y(1,3);error stop 'dde accuracy';end if
  s=ode(rhs,[1.0_dp],[0.0_dp,1.0_dp],method='ode45',rtol=1e-9_dp,atol=1e-11_dp)
  if(.not.s%ok().or.abs(s%y(1,2)-exp(-1.0_dp))>2e-7_dp)error stop 'facade ode45'
  s=ode(rhs,[1.0_dp],[0.0_dp,1.0_dp])
  if(.not.s%ok().or.abs(s%y(1,2)-exp(-1.0_dp))>2e-7_dp)error stop 'facade lsoda'
  print *,'test_dde_facade: PASS'
contains
  subroutine dr(t,y,h,dy)
    use desolve_utilities, only : history_buffer
    real(dp),intent(in)::t,y(:);type(history_buffer),intent(in)::h;real(dp),intent(out)::dy(:)
    real(dp),allocatable::lag(:)
    lag=h%lag_value(t-1.0_dp);dy(1)=-lag(1);if(y(1)<-huge(1.0_dp))stop
  end subroutine dr
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:);real(dp),intent(out)::dy(:);dy=-y;if(t< -huge(1.0_dp))stop
  end subroutine rhs
end program test_dde_facade
