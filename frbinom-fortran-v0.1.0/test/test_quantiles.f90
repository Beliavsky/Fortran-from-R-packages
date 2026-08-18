program test_quantiles
   use frbinom
   implicit none
   real(dp), allocatable :: cdf(:)
   real(dp) :: p
   integer :: q,k,fails,status

   fails = 0
   call frbinom_cdf_table(40,0.55_dp,0.72_dp,0.16_dp,.false.,cdf,status)
   if (status /= 0) fails=fails+1

   do k = 1, 99
      p = real(k,dp)/100.0_dp
      q = qfrbinom(p,40,0.55_dp,0.72_dp,0.16_dp)
      if (cdf(q) < p-5.0e-15_dp) fails=fails+1
      if (q > 0) then
         if (cdf(q-1) >= p+5.0e-15_dp) fails=fails+1
      end if
   end do

   call frbinom2_cdf_table(40,0.8_dp,0.2_dp,0.1_dp,.false.,cdf,status)
   do k = 1, 99
      p = real(k,dp)/100.0_dp
      q = qfrbinom2(p,40,0.8_dp,0.2_dp,0.1_dp)
      if (cdf(q) < p-5.0e-15_dp) fails=fails+1
      if (q > 0) then
         if (cdf(q-1) >= p+5.0e-15_dp) fails=fails+1
      end if
   end do

   if (qfrbinom(0.0_dp,40,0.55_dp,0.72_dp,0.16_dp) /= 0) fails=fails+1
   if (qfrbinom(1.0_dp,40,0.55_dp,0.72_dp,0.16_dp) /= 40) fails=fails+1

   if (fails /= 0) then
      print *, "test_quantiles: FAIL", fails
      error stop 1
   end if
   print *, "test_quantiles: PASS"
end program test_quantiles
