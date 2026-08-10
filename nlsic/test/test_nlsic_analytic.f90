module analytic_problem
   use nlsic, only : dp
   implicit none
contains
   subroutine exp_residual(par,r,ierr)
      real(dp),intent(in)::par(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
      integer::i; real(dp)::x,s
      do i=1,size(r)
         x=real(i-1,dp); s=exp(par(1)*x+par(2)); r(i)=s-exp(x+2.0_dp)
      end do
      ierr=0
   end subroutine
   subroutine exp_jacobian(par,r,j,ierr)
      real(dp),intent(in)::par(:); real(dp),intent(out)::r(:),j(:,:); integer,intent(out)::ierr
      integer::i; real(dp)::x,s
      do i=1,size(r)
         x=real(i-1,dp); s=exp(par(1)*x+par(2)); r(i)=s-exp(x+2.0_dp)
         j(i,1)=s*x; j(i,2)=s
      end do
      ierr=0
   end subroutine
end module analytic_problem

program test_nlsic_analytic
   use nlsic
   use analytic_problem
   implicit none
   type(nlsic_result)::r
   type(nlsic_control)::c
   real(dp)::p0(2),u(2,2),co(2)
   p0=0.0_dp; u=0.0_dp; u(1,1)=1.0_dp; u(2,2)=1.0_dp; co=1.0_dp
   c%history=.true.; c%report_ci=.true.
   call nlsic_solve(p0,6,exp_residual,r,jacobian=exp_jacobian,u=u,co=co,control=c)
   call assert_true(r%succeeded(),'nlsic analytic status')
   call assert_true(r%converged,'nlsic analytic converged')
   call assert_vec(r%par,[1.0_dp,2.0_dp],2.0e-8_dp,'nlsic analytic parameters')
   call assert_true(sum(r%residuals*r%residuals)<1.0e-18_dp,'nlsic analytic residual')
   call assert_true(minval(matmul(u,r%par)-co)>=-1.0e-10_dp,'nlsic analytic constraints')
   call assert_true(size(r%par_history,2)>=2,'history retained')
   call assert_true(all(r%hci>=0.0_dp),'confidence interval widths')
   print *, 'test_nlsic_analytic: PASS'
contains
   subroutine assert_vec(x,y,tol,msg)
      real(dp),intent(in)::x(:),y(:),tol; character(*),intent(in)::msg
      if(maxval(abs(x-y))>tol) then; print *, 'FAIL: ',msg,x,y; error stop 1; end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok) then; print *, 'FAIL: ',msg; error stop 1; end if
   end subroutine
end program test_nlsic_analytic
