module constrained_exponential_problem
   use nlsic, only : dp
   implicit none
contains
   subroutine residual(par,r,ierr)
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: r(:)
      integer, intent(out) :: ierr
      integer :: i
      real(dp) :: x,s
      do i=1,size(r)
         x=real(i-1,dp)
         s=exp(par(1)*x+par(2))
         r(i)=s-exp(x+2.0_dp)
      end do
      ierr=0
   end subroutine residual

   subroutine jacobian(par,r,j,ierr)
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: r(:),j(:,:)
      integer, intent(out) :: ierr
      integer :: i
      real(dp) :: x,s
      do i=1,size(r)
         x=real(i-1,dp)
         s=exp(par(1)*x+par(2))
         r(i)=s-exp(x+2.0_dp)
         j(i,1)=s*x
         j(i,2)=s
      end do
      ierr=0
   end subroutine jacobian
end module constrained_exponential_problem

program constrained_exponential
   use nlsic
   use constrained_exponential_problem
   implicit none
   real(dp) :: par0(2),u(2,2),co(2)
   type(nlsic_result) :: fit
   type(nlsic_control) :: control

   par0=0.0_dp
   u=0.0_dp; u(1,1)=1.0_dp; u(2,2)=1.0_dp
   co=1.0_dp
   control%history=.true.

   call nlsic_solve(par0,6,residual,fit,jacobian=jacobian,u=u,co=co,control=control)
   print '(a,i0)', 'status = ',fit%status
   print '(a,2f16.10)', 'par    = ',fit%par
   print '(a,es16.8)', 'RSS    = ',sum(fit%residuals*fit%residuals)
end program constrained_exponential
