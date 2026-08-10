module bad_problem
   use nlsic, only : dp
   implicit none
contains
   subroutine residual(par,r,ierr)
      real(dp),intent(in)::par(:); real(dp),intent(out)::r(:); integer,intent(out)::ierr
      r=par(1)-1.0_dp; ierr=0
   end subroutine
end module bad_problem

program test_regularized_errors
   use nlsic
   use bad_problem
   implicit none
   type(lsi_result)::lr
   type(nlsic_result)::nr
   type(nlsic_control)::c
   real(dp)::a(4,3),b(4),u(2,2),co(2),p0(2)
   integer::i
   do i=1,4
      a(i,1)=real(i,dp); a(i,2)=1.0_dp; a(i,3)=a(i,2)
   end do
   b=matmul(a,[0.5_dp,1.0_dp,-0.25_dp])
   call lsi(a,b,lr)
   call assert_true(lr%status==LSI_RANK_DEFICIENT,'lsi rank-deficient status')
   call lsi_reg(a,b,lr)
   call assert_true(lr%succeeded(),'lsi_reg status')
   call assert_true(lr%lambda>0.0_dp,'lsi_reg lambda')
   call assert_true(vecnorm(matmul(a,lr%x)-b)<1.0e-5_dp,'lsi_reg residual')

   p0=0.0_dp; u=0.0_dp; u(1,1)=1.0_dp; u(2,1)=-1.0_dp; co=[1.0_dp,1.0_dp]
   c%report_ci=.false.
   call nlsic_solve(p0,2,residual,nr,u=u,co=co,control=c)
   call assert_true(nr%status==NLSIC_INFEASIBLE,'infeasible nonlinear constraints')
   print *, 'test_regularized_errors: PASS'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok) then; print *, 'FAIL: ',msg; error stop 1; end if
   end subroutine
end program test_regularized_errors
