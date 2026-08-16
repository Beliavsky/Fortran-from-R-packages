program test_basic
   use nleqslv_fortran
   implicit none
   type(nleq_result) :: r
   type(nleq_options) :: opt
   real(dp) :: x0(2)
   integer :: fails
   fails=0
   x0=[2.0_dp,0.5_dp]
   opt=nleq_options()
   opt%return_jacobian=.true.
   call solve_nleqslv(x0, f, r, opt)
   if(r%termcd /= 1) fails=fails+1
   if(maxval(abs(r%x-[1.0_dp,1.0_dp])) > 1.0e-6_dp) fails=fails+1
   if(maxval(abs(r%fvec)) > 1.0e-7_dp) fails=fails+1
   opt%method=NLEQ_NEWTON
   call solve_nleqslv(x0, f, r, opt, jac)
   if(r%termcd /= 1) fails=fails+1
   if(maxval(abs(r%x-[1.0_dp,1.0_dp])) > 1.0e-8_dp) fails=fails+1
   if(maxval(abs(r%jac-jacval(r%x))) > 1.0e-10_dp) fails=fails+1
   if(fails==0) then
      print *, 'test_basic: PASS'
   else
      print *, 'test_basic: FAIL',fails
      error stop 1
   end if
contains
   subroutine f(x,y)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      y(1)=x(1)**2+x(2)**2-2.0_dp
      y(2)=exp(x(1)-1.0_dp)+x(2)**3-2.0_dp
   end subroutine
   subroutine jac(x,j)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::j(:,:)
      j=jacval(x)
   end subroutine
   function jacval(x) result(j)
      real(dp),intent(in)::x(:)
      real(dp)::j(2,2)
      j(1,1)=2*x(1); j(1,2)=2*x(2)
      j(2,1)=exp(x(1)-1); j(2,2)=3*x(2)**2
   end function
end program
