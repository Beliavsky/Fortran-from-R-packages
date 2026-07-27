! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
program vrtest_demo
   use vrtest
   implicit none
   integer, parameter :: n = 120
   integer, parameter :: kvec(3) = [2,5,10]
   real(dp) :: y(n)
   integer :: i
   type(auto_vr_result) :: av
   type(auto_q_result) :: aq
   type(lmcd_result) :: lm
   type(chow_denning_result) :: cd
   type(wright_result) :: wr
   type(spectral_shape_result) :: ss

   do i = 1, n
      y(i) = 0.012_dp*sin(0.47_dp*real(i,dp)) + &
         0.007_dp*cos(0.19_dp*real(i,dp)) + &
         0.003_dp*sin(1.13_dp*real(i,dp))
   end do

   av = automatic_variance_ratio(y)
   aq = automatic_portmanteau(y,10)
   lm = lo_mackinlay(y,kvec)
   cd = chow_denning(y,kvec)
   wr = wright_tests(y,kvec)
   ss = spectral_shape_test(y)

   print '(a)', 'vrtest-fortran demonstration'
   print '(a,i0)', 'observations: ',n
   print '(a,f12.6)', 'automatic VR statistic: ',av%statistic
   print '(a,f12.6)', 'automatic VR bandwidth: ',av%bandwidth
   print '(a,f12.6)', 'automatic Q statistic:  ',aq%statistic
   print '(a,i0)', 'automatic Q selected lag: ',aq%selected_lag
   print '(a,f10.6)', 'automatic Q p-value: ',aq%p_value
   print '(a)', 'Lo-MacKinlay statistics (k, M1, M2):'
   do i = 1, size(kvec)
      print '(i4,2f14.6)',kvec(i),lm%homoskedastic(i),lm%heteroskedastic(i)
   end do
   print '(a,2f14.6)', 'Chow-Denning maxima: ',cd%cd_homoskedastic,cd%cd_heteroskedastic
   print '(a)', 'Wright statistics (k, R1, R2, S1):'
   do i = 1, size(kvec)
      print '(i4,3f14.6)',kvec(i),wr%statistics(i,:)
   end do
   print '(a,3f14.6)', 'Spectral shape AD, CvM, M: ', &
      ss%anderson_darling,ss%cramer_von_mises,ss%maximum
end program vrtest_demo
