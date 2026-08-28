module sweep_problem
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
program test_sweep
   use nleqslv_fortran
   use sweep_problem
   implicit none
   type(nleq_test_result)::r
   real(dp)::x0(2)
   integer::fails
   x0=[0.5_dp,0.5_dp]; fails=0
   call test_nleqslv(x0,dsln,r)
   if(size(r%termcd)/=14) fails=fails+1
   if(any(r%termcd/=1 .and. r%termcd/=2)) fails=fails+1
   if(maxval(r%fnorm)>1.0e-12_dp) fails=fails+1
   if(fails==0) then
      print *, 'test_sweep: PASS'
   else
      print *, 'test_sweep: FAIL',fails; error stop 1
   end if
end program
