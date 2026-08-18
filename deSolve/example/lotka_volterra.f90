module lotka_model
  use desolve, only : dp
  implicit none
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::dy(:)
    real(dp),parameter::a=1.5_dp,b=1.0_dp,c=3.0_dp,d=1.0_dp
    dy(1)=a*y(1)-b*y(1)*y(2)
    dy(2)=d*y(1)*y(2)-c*y(2)
    if(t < -huge(1.0_dp))stop
  end subroutine rhs
end module lotka_model

program lotka_volterra
  use desolve, only : dp,ode_result,ode
  use lotka_model, only : rhs
  implicit none
  type(ode_result)::sol
  real(dp)::times(11)
  integer::i
  do i=1,size(times);times(i)=0.5_dp*real(i-1,dp);end do
  sol=ode(rhs,[10.0_dp,5.0_dp],times,method='lsoda',rtol=1e-9_dp,atol=1e-11_dp)
  if(.not.sol%ok())error stop sol%message
  write(*,'(a)')' time       prey       predator'
  do i=1,size(times);write(*,'(f6.2,2(1x,f11.6))')sol%t(i),sol%y(1,i),sol%y(2,i);end do
end program lotka_volterra
