! SPDX-License-Identifier: GPL-2.0-or-later
program test_stode
  use rootsolve, only : dp, stode, steady_band, steady_result, steady_options
  implicit none
  type(steady_result)::a,b,c
  type(steady_options)::opt
  real(dp)::y0(3)
  y0=[4.0_dp,-2.0_dp,7.0_dp]
  opt%rtol=1e-10_dp;opt%atol=1e-12_dp
  a=stode(rhs,y0,options=opt)
  if(.not.a%steady.or.maxval(abs(a%y-1.0_dp))>1e-9_dp)error stop 1
  opt%jactype='fullusr'
  b=stode(rhs,y0,options=opt,jacfunc=jac)
  if(.not.b%steady.or.maxval(abs(b%y-1.0_dp))>1e-10_dp)error stop 2
  c=steady_band(rhs,y0,1,1,options=opt)
  if(.not.c%steady.or.maxval(abs(c%y-1.0_dp))>1e-9_dp)error stop 3
  print *, 'test_stode: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::dy(:)
    dy(1)=2.0_dp*y(1)-y(2)-1.0_dp
    dy(2)=-y(1)+2.0_dp*y(2)-y(3)
    dy(3)=-y(2)+2.0_dp*y(3)-1.0_dp
    if(t < -huge(1.0_dp))error stop 99
  end subroutine rhs
  subroutine jac(t,y,j)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::j(:,:)
    j=0.0_dp
    j(1,1)=2.0_dp;j(1,2)=-1.0_dp
    j(2,1)=-1.0_dp;j(2,2)=2.0_dp;j(2,3)=-1.0_dp
    j(3,2)=-1.0_dp;j(3,3)=2.0_dp
    if(t+sum(y) < -huge(1.0_dp))error stop 99
  end subroutine jac
end program test_stode
