module numeric_problem
   use nlsic, only : dp
   implicit none
contains
   subroutine residual(par,r,ierr)
      real(dp),intent(in)::par(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
      r(1)=par(1)**2-1.0_dp
      r(2)=par(2)-2.0_dp
      r(3)=0.5_dp*(par(1)+par(2)-3.0_dp)
      ierr=0
   end subroutine
end module numeric_problem

program test_nlsic_numeric_equality
   use nlsic
   use numeric_problem
   implicit none
   type(nlsic_result)::r
   type(nlsic_control)::c
   real(dp)::p0(2),e(1,2),eco(1)
   p0=0.0_dp; e=reshape([1.0_dp,1.0_dp],[1,2]); eco=3.0_dp
   c%report_ci=.true.; c%least_norm_step=.true.
   call nlsic_solve(p0,3,residual,r,e=e,eco=eco,control=c)
   call assert_true(r%succeeded(),'numeric Jacobian status')
   call assert_vec(r%par,[1.0_dp,2.0_dp],5.0e-7_dp,'numeric Jacobian/equality parameters')
   call assert_close(dot_product(e(1,:),r%par),3.0_dp,1.0e-9_dp,'equality exact')
   call assert_true(sum(r%residuals*r%residuals)<1.0e-12_dp,'numeric residual')
   print *, 'test_nlsic_numeric_equality: PASS'
contains
   subroutine assert_close(x,y,tol,msg)
      real(dp),intent(in)::x,y,tol; character(*),intent(in)::msg
      if(abs(x-y)>tol) then; print *, 'FAIL: ',msg,x,y; error stop 1; end if
   end subroutine
   subroutine assert_vec(x,y,tol,msg)
      real(dp),intent(in)::x(:),y(:),tol; character(*),intent(in)::msg
      if(maxval(abs(x-y))>tol) then; print *, 'FAIL: ',msg,x,y; error stop 1; end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok) then; print *, 'FAIL: ',msg; error stop 1; end if
   end subroutine
end program test_nlsic_numeric_equality
