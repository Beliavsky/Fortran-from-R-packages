! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
program basic_tests
   use vrtest
   implicit none
   real(dp), parameter :: returns(16) = [ &
      0.011_dp,-0.004_dp,0.008_dp,0.002_dp,-0.010_dp,0.006_dp,0.014_dp,-0.007_dp, &
      0.003_dp,0.009_dp,-0.005_dp,0.004_dp,-0.012_dp,0.007_dp,0.005_dp,-0.002_dp ]
   integer, parameter :: periods(3) = [2,4,8]
   type(lmcd_result) :: result
   integer :: i

   result = lo_mackinlay(returns,periods)
   print '(a)', 'holding period       M1            M2'
   do i = 1, size(periods)
      print '(i8,2f14.6)',periods(i),result%homoskedastic(i),result%heteroskedastic(i)
   end do
end program basic_tests
