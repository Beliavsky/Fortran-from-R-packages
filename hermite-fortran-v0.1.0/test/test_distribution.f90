program test_distribution
   use hermite
   implicit none
   real(dp) :: p,s
   integer :: k,fails

   fails=0
   if (abs(dhermite_exact(0.0_dp,0.8_dp,0.3_dp,3)- &
       0.33287108369807955_dp) > 2.0e-14_dp) fails=fails+1
   if (abs(dhermite_exact(3.0_dp,0.8_dp,0.3_dp,3)- &
       0.12826632425165999_dp) > 2.0e-14_dp) fails=fails+1
   if (abs(dhermite_exact(7.0_dp,0.8_dp,0.3_dp,3)- &
       0.01370150978029439_dp) > 2.0e-14_dp) fails=fails+1

   if (abs(int_hermite(3_i64,0.8_dp,0.3_dp,3)- &
       dhermite_exact(3.0_dp,0.8_dp,0.3_dp,3)) > 2.0e-14_dp) fails=fails+1

   s=0.0_dp
   do k=0,60
      s=s+dhermite_exact(real(k,dp),0.8_dp,0.3_dp,3)
   end do
   if (abs(s-1.0_dp) > 3.0e-14_dp) fails=fails+1

   p=phermite(4.0_dp,0.8_dp,0.3_dp,3,exact=.true.)
   if (abs(p-0.919523081607575_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(hermite_mean(0.8_dp,0.3_dp,3)-1.7_dp) > 1.0e-14_dp) fails=fails+1
   if (abs(hermite_variance(0.8_dp,0.3_dp,3)-3.5_dp) > 1.0e-14_dp) fails=fails+1

   if (fails/=0) then
      print *,'test_distribution: FAIL',fails
      error stop 1
   end if
   print *,'test_distribution: PASS'
end program test_distribution
