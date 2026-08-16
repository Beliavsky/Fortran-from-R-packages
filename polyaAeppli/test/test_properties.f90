program test_properties
   use polya_aeppli
   implicit none
   integer :: fails,k,q
   real(dp) :: s,cprev,c,p
   real(dp), allocatable :: lp(:),lc(:),ls(:)
   integer :: status

   fails = 0

   s = 0.0_dp
   do k = 0, 200
      s = s+d_polya_aeppli(real(k,dp),3.2_dp,0.4_dp)
   end do
   if (abs(s-1.0_dp) > 2.0e-13_dp) fails=fails+1

   cprev = 0.0_dp
   do k = 0, 50
      c = p_polya_aeppli(real(k,dp),3.2_dp,0.4_dp)
      if (c < cprev-2.0e-15_dp) fails=fails+1
      if (abs((c-cprev)-d_polya_aeppli(real(k,dp),3.2_dp,0.4_dp)) > 3.0e-14_dp) &
         fails=fails+1
      cprev = c
   end do

   do k = 1, 99
      p = real(k,dp)/100.0_dp
      q = q_polya_aeppli(p,3.2_dp,0.4_dp)
      if (p_polya_aeppli(real(q,dp),3.2_dp,0.4_dp) < p-5.0e-15_dp) fails=fails+1
      if (q > 0) then
         if (p_polya_aeppli(real(q-1,dp),3.2_dp,0.4_dp) >= p+5.0e-15_dp) &
            fails=fails+1
      end if
   end do

   call log_pmf_array(80,3.2_dp,0.4_dp,lp,status)
   call log_cdf_array(lp,lc)
   call log_sf_array(80,3.2_dp,0.4_dp,lp,ls,status)
   do k = 0, 80
      if (abs(lc(k)-p_polya_aeppli(real(k,dp),3.2_dp,0.4_dp,.true.,.true.)) > &
          2.0e-13_dp) fails=fails+1
      if (abs(ls(k)-p_polya_aeppli(real(k,dp),3.2_dp,0.4_dp,.false.,.true.)) > &
          2.0e-12_dp) fails=fails+1
   end do

   if (abs(d_polya_aeppli(1.5_dp,3.2_dp,0.4_dp)) > epsilon(1.0_dp)) fails=fails+1

   if (fails /= 0) then
      print *, "test_properties: FAIL", fails
      error stop 1
   end if
   print *, "test_properties: PASS"
end program test_properties
