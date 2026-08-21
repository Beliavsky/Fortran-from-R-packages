! SPDX-License-Identifier: GPL-2.0-or-later
program test_pde_wrappers
  use rootsolve, only : dp, steady_2d, steady_3d, steady_result, steady_options
  implicit none
  type(steady_result)::a,b
  type(steady_options)::op
  real(dp)::y2(24),y3(24)
  y2=3.0_dp;y3=-2.0_dp
  op%rtol=1e-10_dp;op%atol=1e-12_dp
  a=steady_2d(diag_rhs,y2,2,[4,3],options=op)
  if(.not.a%steady.or.maxval(abs(a%y-1.0_dp))>1e-9_dp)error stop 1
  b=steady_3d(diag_rhs,y3,1,[4,3,2],options=op)
  if(.not.b%steady.or.maxval(abs(b%y-1.0_dp))>1e-9_dp)error stop 2
  print *, 'test_pde_wrappers: PASS'
contains
  subroutine diag_rhs(t,y,dy)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::dy(:)
    dy=y-1.0_dp
    if(t < -huge(1.0_dp))error stop 99
  end subroutine diag_rhs
end program test_pde_wrappers
