program basic_usage
   use grpreg
   implicit none
   integer, parameter :: n = 12, p = 4
   real(dp) :: x(n,p), y(n)
   integer :: group(p), i
   type(grpreg_fit_type) :: fit
   real(dp), allocatable :: pred(:)
   do i = 1, n
      x(i,1) = real(i,dp)/real(n,dp)
      x(i,2) = sin(real(i,dp))
      x(i,3) = cos(0.7_dp*real(i,dp))
      x(i,4) = real(mod(i,3),dp)
   end do
   group = [1,1,2,2]
   y = 1.0_dp+2.0_dp*x(:,1)-x(:,2)+0.5_dp*x(:,3)
   call grpreg_fit(x,y,group,fit,penalty='grMCP',nlambda=20)
   call predict_grpreg(fit,x,fit%lambda(fit%nlambda),pred,type='response')
   write(*,'(a,i0)') 'Number of lambda values: ',fit%nlambda
   write(*,'(a,es12.4)') 'Smallest-lambda RSS: ',sum((y-pred)**2)
end program basic_usage
