! SPDX-License-Identifier: GPL-3.0-only
program derivative_example
   use smoots
   implicit none
   integer,parameter::n=180
   real(dp)::y(n),x
   integer::i
   type(smooth_result)::d1,d2
   do i=1,n
      x=real(i,dp)/real(n,dp)
      y(i)=exp(x)+0.02_dp*sin(31.0_dp*x)
   end do
   call dsmooth(y,d1,d=1)
   call dsmooth(y,d2,d=2)
   print '(a,f10.6)', 'first derivative bandwidth: ',d1%b0
   print '(a,f10.6)', 'second derivative bandwidth:',d2%b0
   print '(a,2f12.6)', 'middle estimates: ',d1%estimate(n/2),d2%estimate(n/2)
end program derivative_example
