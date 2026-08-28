module strategy_problem
   use nleqslv_fortran, only : dp
   implicit none
contains
   subroutine dsln(x,y)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      y(1)=x(1)**2+x(2)**2-2.0_dp
      y(2)=exp(x(1)-1.0_dp)+x(2)**3-2.0_dp
   end subroutine
end module
program test_strategies
   use nleqslv_fortran
   use strategy_problem
   implicit none
   type(nleq_options)::o
   type(nleq_result)::r
   real(dp)::x0(2)
   integer::m,g,fails
   fails=0; x0=[0.5_dp,0.5_dp]
   do m=NLEQ_NEWTON,NLEQ_BROYDEN
      do g=NLEQ_NONE,NLEQ_HOOK
         o=nleq_options(); o%method=m; o%global=g
         call solve_nleqslv(x0,dsln,r,o)
         if(r%termcd/=1 .and. r%termcd/=2) then
            print *, 'unexpected termcd',m,g,r%termcd
            fails=fails+1
         end if
         if(maxval(abs(r%fvec))>1.0e-6_dp) fails=fails+1
      end do
   end do
   if(fails==0) then
      print *, 'test_strategies: PASS'
   else
      print *, 'test_strategies: FAIL',fails; error stop 1
   end if
end program
