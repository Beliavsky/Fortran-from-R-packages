program basic_vgam
   use vgam
   implicit none
   integer, parameter :: n=12
   real(dp) :: x(n,2), y(n), z(n)
   type(vglm_result_t) :: fit
   type(vgam_smooth_result_t) :: smooth_fit
   integer :: i

   do i=1,n
      z(i)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      x(i,:)=[1.0_dp,z(i)]
      y(i)=exp(0.25_dp+0.60_dp*z(i))
   end do

   call fit_poisson(y,x,fit)
   print '(a,2f12.6)', 'Poisson coefficients: ',fit%coefficients

   y=sin(pi*z)+0.05_dp*cos(6.0_dp*z)
   call fit_gam_gaussian(z,y,smooth_fit,df=6,lambda=0.2_dp)
   print '(a,l1)', 'Spline fit converged: ',smooth_fit%fit%converged
   print '(a,f12.6)', 'Spline residual deviance: ',smooth_fit%fit%deviance
end program basic_vgam
