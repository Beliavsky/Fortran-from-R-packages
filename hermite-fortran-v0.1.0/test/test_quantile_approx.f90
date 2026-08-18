program test_quantile_approx
   use hermite
   implicit none
   integer(i64) :: q
   real(dp) :: p,c,e
   integer :: k,fails

   fails=0
   do k=1,99
      p=real(k,dp)/100.0_dp
      q=qhermite(p,2.5_dp,1.2_dp,3,exact=.true.)
      c=phermite(real(q,dp),2.5_dp,1.2_dp,3,exact=.true.)
      if (c < p-5.0e-14_dp) fails=fails+1
      if (q>0_i64) then
         if (phermite(real(q-1_i64,dp),2.5_dp,1.2_dp,3,exact=.true.) >= &
             p+5.0e-14_dp) fails=fails+1
      end if
   end do

   p=phermite(120.0_dp,70.0_dp,12.0_dp,3)
   if (p < 0.0_dp .or. p > 1.0_dp) fails=fails+1
   q=qhermite(0.9_dp,70.0_dp,12.0_dp,3)
   if (q < 0_i64) fails=fails+1
   e=edg(120.0_dp,70.0_dp,12.0_dp,3)
   if (abs(p-min(1.0_dp,max(0.0_dp,e))) > 2.0e-14_dp) fails=fails+1

   if (fails/=0) then
      print *,'test_quantile_approx: FAIL',fails
      error stop 1
   end if
   print *,'test_quantile_approx: PASS'
end program test_quantile_approx
