! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
program bootstrap_tests
   use vrtest
   implicit none
   integer, parameter :: n = 100
   real(dp) :: returns(n)
   integer :: i
   type(auto_bootstrap_result) :: result

   do i = 1, n
      returns(i) = 0.01_dp*sin(0.53_dp*real(i,dp)) + &
         0.004_dp*cos(0.17_dp*real(i,dp))
   end do
   call seed_random(12345)
   result = automatic_vr_bootstrap(returns,200,'Mammen')

   print '(a,f12.6)', 'automatic VR statistic: ',result%test_statistic
   print '(a,f10.6)', 'bootstrap p-value: ',result%p_value
   print '(a,2f12.6)', '95 percent statistic interval: ',result%statistic_interval
end program bootstrap_tests
