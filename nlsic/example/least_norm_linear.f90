program least_norm_linear
   use nlsic
   implicit none
   real(dp) :: a(2,3),b(2),u(2,3),co(2)
   type(lsi_result) :: fit

   a=0.0_dp
   a(1,1)=1.0_dp
   a(2,2)=1.0_dp
   b=0.0_dp
   u=0.0_dp
   u(1,:)=1.0_dp
   u(2,3)=-1.0_dp
   co=[1.0_dp,-0.5_dp]

   call lsi_ln(a,b,fit,u,co)
   print '(a,i0)', 'status = ',fit%status
   print '(a,3f14.8)', 'x      = ',fit%x
   print '(a,es14.6)', 'RSS    = ',fit%objective
end program least_norm_linear
