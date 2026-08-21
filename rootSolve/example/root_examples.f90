! SPDX-License-Identifier: GPL-2.0-or-later
program root_examples
  use rootsolve, only : dp, multiroot, root_result, stode, steady_result
  implicit none
  type(root_result)::r
  type(steady_result)::s
  r=multiroot(equations,[1.0_dp,1.0_dp])
  write(*,'(a,2f14.8)') 'multiroot: ',r%root
  s=stode(rhs,[5.0_dp,-2.0_dp])
  write(*,'(a,2f14.8)') 'steady state: ',s%y
contains
  subroutine equations(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    f(1)=x(1)**2+x(2)**2-1.0_dp
    f(2)=x(1)**2-x(2)**2+0.5_dp
  end subroutine equations
  subroutine rhs(t,y,dy)
    real(dp),intent(in)::t,y(:)
    real(dp),intent(out)::dy(:)
    dy(1)=1.0_dp-y(1)
    dy(2)=2.0_dp-2.0_dp*y(2)
    if(t < -huge(1.0_dp))error stop 99
  end subroutine rhs
end program root_examples
