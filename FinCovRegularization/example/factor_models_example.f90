! SPDX-License-Identifier: GPL-2.0-only
program factor_models_example
   use fincovregularization
   implicit none
   real(dp) :: factors(12,2), assets(12,3), estimate(3,3)
   integer :: i, status

   do i = 1, 12
      factors(i,1) = 0.01_dp*real(i-6,dp)
      factors(i,2) = 0.008_dp*real(mod(5*i,11)-5,dp)
   end do
   assets(:,1) = 0.002_dp + 1.1_dp*factors(:,1) - 0.2_dp*factors(:,2)
   assets(:,2) = -0.001_dp + 0.4_dp*factors(:,1) + 0.9_dp*factors(:,2)
   assets(:,3) = 0.003_dp - 0.3_dp*factors(:,1) + 0.7_dp*factors(:,2)

   estimate = macro_factor_cov(assets, factors, status)
   if (status /= fincov_ok) error stop fincov_status_message(status)

   print '(a)', 'Macroeconomic factor covariance estimate:'
   do i = 1, 3
      print '(3es14.5)', estimate(i,:)
   end do
end program factor_models_example
