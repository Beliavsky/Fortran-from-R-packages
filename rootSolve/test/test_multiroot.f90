! SPDX-License-Identifier: GPL-2.0-or-later
program test_multiroot
  use rootsolve, only : dp, multiroot, multiroot_1d, root_result, steady_options
  implicit none
  type(root_result) :: r, r1
  type(steady_options) :: opt
  opt%rtol=1.0e-10_dp
  opt%atol=1.0e-12_dp
  r=multiroot(model,[1.0_dp,1.0_dp],opt)
  if(.not.r%converged) error stop 1
  if(maxval(abs(r%root-[0.5_dp,sqrt(0.75_dp)]))>1e-8_dp) error stop 2
  r1=multiroot_1d(model,[1.0_dp,1.0_dp],nspec=1,options=opt)
  if(.not.r1%converged.or.maxval(abs(r1%f_root))>1e-8_dp) error stop 3
  print *, 'test_multiroot: PASS'
contains
  subroutine model(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    f(1)=x(1)**2+x(2)**2-1.0_dp
    f(2)=x(1)**2-x(2)**2+0.5_dp
  end subroutine model
end program test_multiroot
