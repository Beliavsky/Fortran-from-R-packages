program example_ordering
   use portvine, only : dp, greedy_dvine_order
   implicit none
   real(dp) :: u(4,100),x
   integer :: i,order(4),status
   do i=1,100
      x=real(i,dp)/101.0_dp
      u(1,i)=x
      u(2,i)=max(1.0e-5_dp,min(1.0_dp-1.0e-5_dp,x+0.04_dp*sin(11.0_dp*x)))
      u(3,i)=real(mod(31*i,101)+1,dp)/102.0_dp
      u(4,i)=max(1.0e-5_dp,min(1.0_dp-1.0e-5_dp,u(3,i)+0.03_dp*cos(7.0_dp*x)))
   end do
   call greedy_dvine_order(u,order,[3],3,status)
   print '(a,4(i0,1x))','D-vine order: ',order
end program example_ordering
