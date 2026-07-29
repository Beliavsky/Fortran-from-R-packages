! SPDX-License-Identifier: GPL-2.0-or-later
program rm2006_demo
   use rm2006_kinds, only : dp
   use rm2006_module, only : rm2006_covariance, rm2006_status_message
   implicit none

   real(dp) :: returns(8, 3)
   real(dp), allocatable :: covariance(:, :, :)
   integer :: status
   integer :: i

   returns = reshape([ &
       0.010_dp, -0.012_dp,  0.006_dp, &
      -0.008_dp,  0.004_dp,  0.009_dp, &
       0.015_dp,  0.011_dp, -0.005_dp, &
      -0.004_dp, -0.007_dp,  0.003_dp, &
       0.006_dp,  0.002_dp,  0.008_dp, &
      -0.011_dp,  0.014_dp, -0.006_dp, &
       0.009_dp, -0.003_dp,  0.010_dp, &
       0.002_dp,  0.005_dp, -0.004_dp  &
      ], shape(returns), order=[2, 1])

   call rm2006_covariance(returns, covariance, status=status)
   if (status /= 0) then
      error stop rm2006_status_message(status)
   end if

   print '(a)', 'RM2006 next-period conditional covariance matrix:'
   do i = 1, size(covariance, 1)
      print '(*(es14.6,1x))', covariance(i, :, size(covariance, 3))
   end do
end program rm2006_demo
