! SPDX-License-Identifier: LGPL-3.0-or-later
program optimal_sharpe
   use sharper, only: dp, sropt_result
   use sharper, only: fit_sropt, infer_sropt, sropt_confint, sric
   implicit none
   real(dp) :: returns(100,3), ci(2)
   type(sropt_result) :: z
   integer :: i

   do i = 1, 100
      returns(i,1) = 0.004_dp+0.02_dp*sin(0.21_dp*real(i,dp))
      returns(i,2) = 0.003_dp+0.015_dp*cos(0.17_dp*real(i,dp))
      returns(i,3) = 0.002_dp+0.01_dp*sin(0.37_dp*real(i,dp)) + &
                     0.006_dp*cos(0.13_dp*real(i,dp))
   end do
   z = fit_sropt(returns,ope=12.0_dp)
   ci = sropt_confint(z,0.90_dp)
   print '(a,f10.5)', 'sample optimal Sharpe: ',z%value
   print '(a,f10.5)', 'KRS inferred SNR:      ',infer_sropt(z,'krs')
   print '(a,f10.5)', 'SRIC:                  ',sric(z)
   print '(a,2f10.5)', '90% confidence limits: ',ci
end program optimal_sharpe
