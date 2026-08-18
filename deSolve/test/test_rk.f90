program test_rk
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result
  use desolve_rk, only : rk_method, rk_method_by_name, rk_integrate
  implicit none
  character(len=8), parameter :: names(19)=[character(len=8) :: &
    'euler','rk2','rk4','rk23','rk23bs','rk34f','rk45f','rk45ck','rk45e', &
    'rk45dp6','rk45dp7','rk78dp','rk78f','irk3r','irk5r','irk4hh','irk6kb','irk4l','irk6l']
  type(rk_method)::m
  type(ode_result)::sol
  real(dp)::y0(1),tt(2),tol
  integer::i
  y0=1.0_dp;tt=[0.0_dp,1.0_dp]
  do i=1,size(names)
    m=rk_method_by_name(trim(names(i)))
    if(m%adaptive)then
      sol=rk_integrate(rhs,y0,tt,m,h=0.05_dp,rtol=1e-8_dp,atol=1e-10_dp)
      tol=3e-6_dp
    else if(m%implicit)then
      sol=rk_integrate(rhs,y0,tt,m,h=0.05_dp)
      tol=5e-5_dp
    else
      sol=rk_integrate(rhs,y0,tt,m,h=0.001_dp)
      tol=5e-4_dp
    end if
    if(.not.sol%ok())then
      print *, trim(names(i)), sol%status, sol%message
      error stop 'rk status'
    end if
    if(abs(sol%y(1,2)-exp(-1.0_dp))>tol)then
      print *, trim(names(i)), sol%y(1,2), exp(-1.0_dp), tol
      error stop 'rk accuracy'
    end if
  end do
  print *, 'test_rk: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:);real(dp),intent(out)::dy(:)
    dy=-y
    if(t < -huge(1.0_dp))stop
  end subroutine rhs
end program test_rk
