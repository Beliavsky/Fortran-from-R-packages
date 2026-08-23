! SPDX-License-Identifier: GPL-3.0-or-later
module test_calculus_functions
   use pracma_kinds, only : dp
   implicit none
contains
   function fsin(x) result(y)
      real(dp),intent(in)::x; real(dp)::y; y=sin(x)
   end function fsin
   function froot(x) result(y)
      real(dp),intent(in)::x; real(dp)::y; y=x*x-2.0_dp
   end function froot
   function fscalar(x) result(y)
      real(dp),intent(in)::x; real(dp)::y; y=(x-1.5_dp)**2+2.0_dp
   end function fscalar
   function fvec(x) result(y)
      real(dp),intent(in)::x(:); real(dp)::y
      y=(x(1)-1.0_dp)**2+2.0_dp*(x(2)+2.0_dp)**2
   end function fvec
end module test_calculus_functions

program test_calculus_optimization
   use pracma
   use test_calculus_functions
   implicit none
   type(quadrature_result)::qi
   type(root_result)::rr
   type(optimization_result)::ob,om
   real(dp)::g(2)
   integer::st
   qi=integral(fsin,0.0_dp,pi_dp)
   call check(qi%converged.and.abs(qi%value-2.0_dp)<1e-9_dp,'integral')
   rr=brentDekker(froot,0.0_dp,2.0_dp)
   call check(rr%converged.and.abs(rr%root-sqrt(2.0_dp))<1e-10_dp,'brent')
   ob=fminbnd(fscalar,-3.0_dp,4.0_dp)
   call check(ob%converged.and.abs(ob%x(1)-1.5_dp)<1e-7_dp,'fminbnd')
   om=fminsearch(fvec,[0.0_dp,0.0_dp],max_iter=2000)
   call check(om%converged.and.maxval(abs(om%x-[1.0_dp,-2.0_dp]))<3e-5_dp,'fminsearch')
   call grad(fvec,[1.5_dp,-1.0_dp],g,status=st)
   call check(st==pracma_ok.and.maxval(abs(g-[1.0_dp,4.0_dp]))<1e-5_dp,'gradient')
   print '(a)','test_calculus_optimization: PASS'
contains
   subroutine check(ok,name)
      logical,intent(in)::ok; character(len=*),intent(in)::name
      if(.not.ok)then; write(*,'(a,1x,a)')'FAIL:',trim(name); error stop 1; end if
   end subroutine check
end program test_calculus_optimization
