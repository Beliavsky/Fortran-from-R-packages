program v02_extended
   use gamlss
   implicit none
   integer,parameter::n=80
   real(dp)::x(n),y(n),w(n)
   real(dp),allocatable::fit(:)
   type(fp_spec_t)::fp
   integer::i,status

   w=1.0_dp
   do i=1,n
      x(i)=0.5_dp+3.5_dp*real(i-1,dp)/real(n-1,dp)
      y(i)=1.0_dp+2.0_dp/x(i)+0.02_dp*sin(real(i,dp))
   end do
   call select_fractional_polynomial(x,y,w,1,fp,fit,status)
   print '(a,f7.2)', 'Selected fractional-polynomial power: ',fp%powers(1)
   print '(a,f10.6)', 'Weighted residual deviance: ',fp%deviance
end program v02_extended
