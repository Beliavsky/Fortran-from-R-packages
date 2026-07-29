! SPDX-License-Identifier: GPL-3.0-only
program fixed_smoothers
   use smoots
   implicit none
   real(dp) :: y(21)
   integer :: i
   type(smooth_result) :: lp,kr
   do i=1,21
      y(i)=real(i,dp)**2+0.1_dp*(-1.0_dp)**i
   end do
   call fixed_gsmooth(y,0,1,1,0.2_dp,1,lp)
   call fixed_knsmooth(y,1,0.2_dp,0,kr)
   print '(a,3f12.6)', 'local polynomial: ',lp%estimate(1),lp%estimate(11),lp%estimate(21)
   print '(a,3f12.6)', 'kernel:          ',kr%estimate(1),kr%estimate(11),kr%estimate(21)
end program fixed_smoothers
