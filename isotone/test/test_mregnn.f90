program test_mregnn
   use isotone
   implicit none
   real(dp) :: x(10,3), y(10)
   type(mregnn_result) :: r, rp
   integer :: i
   do i=1,10
      x(i,1)=real(i,dp)
      x(i,2)=real(i*i,dp)
      x(i,3)=real(i*i*i,dp)
   end do
   do i=1,3
      x(:,i)=x(:,i)-sum(x(:,i))/10.0_dp
      x(:,i)=x(:,i)/sqrt(sum(x(:,i)**2))
   end do
   y=x(:,1)+x(:,2)+x(:,3)+[0.1_dp,-0.2_dp,0.05_dp,0.3_dp,-0.1_dp, &
      0.0_dp,-0.15_dp,0.2_dp,-0.05_dp,0.1_dp]
   call mregnn_monotone(x,y,r)
   if(r%status/=0) error stop 'mregnn monotone failed'
   do i=1,9
      if(r%xb(i)>r%xb(i+1)+1.0e-9_dp) error stop 'mregnn monotonicity failed'
   end do
   call mregnn_positive(abs(x)+0.1_dp,abs(y)+0.2_dp,rp)
   if(rp%status/=0) error stop 'mregnn positive failed'
   if(minval(rp%xb)<-1.0e-9_dp) error stop 'mregnn positivity failed'
   print *, 'test_mregnn: PASS'
end program
