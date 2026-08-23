program test_random
   use psdistr, only : dp, rtppn, rpc, rdsn, ren, rspc, reck
   implicit none
   integer, parameter :: n=20000
   real(dp) :: x(n), m
   integer :: fails
   fails=0

   call rtppn(n,0.0_dp,1.0_dp,1.0_dp,1.0_dp,x)
   m=sum(x)/real(n,dp)
   if(abs(m)>0.04_dp) fails=fails+1

   call rpc(n,0.0_dp,1.0_dp,1.0_dp,x)
   m=sum(x)/real(n,dp)
   if(abs(m)>0.04_dp) fails=fails+1

   call rdsn(n,0.0_dp,1.0_dp,0.0_dp,0.0_dp,x)
   m=sum(x)/real(n,dp)
   if(abs(m)>0.04_dp) fails=fails+1

   call ren(n,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,x)
   if(any(.not.(x < huge(1.0_dp)))) fails=fails+1

   call rspc(n,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,x)
   m=sum(x)/real(n,dp)
   if(abs(m)>0.04_dp) fails=fails+1

   call reck(n,2.0_dp,1.0_dp,x)
   if(any(abs(x)>2.0_dp+1e-12_dp)) fails=fails+1
   m=sum(x)/real(n,dp)
   if(abs(m)>0.04_dp) fails=fails+1

   if(fails/=0) then
      write(*,'(a,i0)') 'test_random: FAIL ',fails
      error stop 1
   end if
   write(*,'(a)') 'test_random: PASS'
end program test_random
