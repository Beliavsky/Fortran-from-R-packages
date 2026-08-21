! SPDX-License-Identifier: GPL-2.0-or-later
program test_derivatives
  use rootsolve, only : dp, gradient, hessian, jacobian_full, jacobian_band
  implicit none
  real(dp) :: x(2), g(2,2), h(2,2), jf(3,3), jb(3,3), y(3)
  x=[1.0_dp,2.0_dp]
  call gradient(vfun,x,g,centered=.true.,pert=1.0e-6_dp)
  if(maxval(abs(g-reshape([2.0_dp,1.0_dp,1.0_dp,4.0_dp],[2,2])))>2.0e-6_dp) error stop 1
  x=[1.0_dp,1.0_dp]
  call hessian(banana,x,h,centered=.true.,pert=1.0e-4_dp)
  if(maxval(abs(h-reshape([802.0_dp,-400.0_dp,-400.0_dp,200.0_dp],[2,2])))>0.2_dp) error stop 2
  y=[1.0_dp,2.0_dp,3.0_dp]
  call jacobian_full(y,rhs,jf,pert=1.0e-7_dp)
  call jacobian_band(y,rhs,1,1,jb,pert=1.0e-7_dp)
  if(abs(jf(1,1)+2.0_dp)>1e-6_dp.or.abs(jf(1,2)-1.0_dp)>1e-6_dp) error stop 3
  if(abs(jb(2,2)+2.0_dp)>1e-6_dp.or.abs(jb(1,2)-1.0_dp)>1e-6_dp) error stop 4
  print *, 'test_derivatives: PASS'
contains
  subroutine vfun(z,f)
    real(dp),intent(in)::z(:)
    real(dp),intent(out)::f(:)
    f(1)=z(1)**2+z(2)
    f(2)=z(1)+z(2)**2
  end subroutine vfun
  real(dp) function banana(z) result(v)
    real(dp),intent(in)::z(:)
    v=100.0_dp*(z(2)-z(1)**2)**2+(1.0_dp-z(1))**2
  end function banana
  subroutine rhs(t,z,dz)
    real(dp),intent(in)::t,z(:)
    real(dp),intent(out)::dz(:)
    dz(1)=-2*z(1)+z(2)
    dz(2)=z(1)-2*z(2)+z(3)
    dz(3)=z(2)-2*z(3)
    if(t < -huge(1.0_dp)) error stop 99
  end subroutine rhs
end program test_derivatives
