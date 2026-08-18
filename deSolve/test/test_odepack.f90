program test_odepack
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result
  use desolve_odepack, only : lsoda,lsode,vode
  implicit none
  type(ode_result)::a,b,c
  real(dp)::tt(3)
  tt=[0.0_dp,0.5_dp,1.0_dp]
  a=lsoda(rhs,[1.0_dp],tt,rtol=1e-10_dp,atol=1e-12_dp)
  b=lsode(rhs,[1.0_dp],tt,rtol=1e-10_dp,atol=1e-12_dp,mf=10)
  c=vode(rhs,[1.0_dp],tt,rtol=1e-10_dp,atol=1e-12_dp,mf=10)
  if(.not.a%ok().or..not.b%ok().or..not.c%ok())error stop 'odepack family status'
  if(max(abs(a%y(1,3)-exp(-1.0_dp)),abs(b%y(1,3)-exp(-1.0_dp)),abs(c%y(1,3)-exp(-1.0_dp)))>3e-8_dp)then
    print *,a%y(1,3),b%y(1,3),c%y(1,3);error stop 'odepack family accuracy'
  end if
  print *,'test_odepack: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:);real(dp),intent(out)::dy(:);dy=-y;if(t< -huge(1.0_dp))stop
  end subroutine rhs
end program test_odepack
