! SPDX-License-Identifier: GPL-2.0-or-later
program custom_parameters
   use rm2006_kinds, only : dp
   use rm2006_module, only : rm2006_covariance
   implicit none

   real(dp) :: returns(5, 2)
   real(dp), allocatable :: covariance(:, :, :)
   integer :: t

   returns = reshape([ &
       0.01_dp, -0.02_dp, &
       0.03_dp,  0.01_dp, &
      -0.02_dp,  0.04_dp, &
       0.00_dp, -0.01_dp, &
       0.05_dp,  0.02_dp  &
      ], shape(returns), order=[2, 1])

   call rm2006_covariance(returns, covariance, tau0=20.0_dp, &
                          tau1=2.0_dp, kmax=3, rho=1.5_dp)

   do t = 1, size(covariance, 3)
      print '(a,i0)', 'forecast index ', t
      print '(*(es14.6,1x))', covariance(:, :, t)
   end do
end program custom_parameters
