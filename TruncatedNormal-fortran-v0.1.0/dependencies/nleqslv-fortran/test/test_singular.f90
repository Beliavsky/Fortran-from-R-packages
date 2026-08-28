module singular_problem
   use nleqslv_fortran, only : dp
   implicit none
contains
   subroutine f(x,y)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
      y(1)=x(1)+x(2)-x(1)*x(2)-2.0_dp
      y(2)=x(1)+x(3)-x(1)*x(3)-3.0_dp
      y(3)=x(2)+x(3)-4.0_dp
   end subroutine
   subroutine jac(x,j)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::j(:,:)
      j(:,1)=[1.0_dp-x(2),1.0_dp-x(3),0.0_dp]
      j(:,2)=[1.0_dp-x(1),0.0_dp,1.0_dp]
      j(:,3)=[0.0_dp,1.0_dp-x(1),1.0_dp]
   end subroutine
end module
program test_singular
   use nleqslv_fortran
   use singular_problem
   implicit none
   type(nleq_options)::o
   type(nleq_result)::r
   real(dp)::x0(3),want(3)
   integer::fails
   fails=0; x0=[1.0_dp,2.0_dp,3.0_dp]; want=[-0.5_dp,5.0_dp/3.0_dp,7.0_dp/3.0_dp]
   o=nleq_options(); o%method=NLEQ_NEWTON; o%allow_singular=.true.
   call solve_nleqslv(x0,f,r,o,jac)
   if(r%termcd/=1) fails=fails+1
   if(maxval(abs(r%x-want))>1.0e-7_dp) fails=fails+1
   if(fails==0) then
      print *, 'test_singular: PASS'
   else
      print *, 'test_singular: FAIL',fails,r%termcd,r%x; error stop 1
   end if
end program
