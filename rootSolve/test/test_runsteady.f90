! SPDX-License-Identifier: GPL-2.0-or-later
program test_runsteady
  use rootsolve, only : dp, runsteady, runsteady_options, steady, steady_options, steady_result
  implicit none
  type(runsteady_options)::op
  type(steady_result)::r,s
  type(steady_options)::sop
  op%stol=1e-9_dp;op%rtol=1e-9_dp;op%atol=1e-11_dp;op%maxsteps=10000
  r=runsteady(rhs,[0.0_dp],[0.0_dp,100.0_dp],op)
  if(.not.r%steady)error stop 1
  if(abs(r%y(1)-3.0_dp)>1e-7_dp)error stop 2
  if(r%estimated_precision>=op%stol)error stop 3
  sop%rtol=1e-9_dp;sop%atol=1e-11_dp
  s=steady(rhs,[0.0_dp],time=0.0_dp,method='runsteady',options=sop)
  if(.not.s%steady.or.abs(s%y(1)-3.0_dp)>1e-7_dp)error stop 4
  print *, 'test_runsteady: PASS'
contains
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::dy(:)
    dy(1)=-0.5_dp*(y(1)-3.0_dp)
    if(t < -huge(1.0_dp))error stop 99
  end subroutine rhs
end program test_runsteady
