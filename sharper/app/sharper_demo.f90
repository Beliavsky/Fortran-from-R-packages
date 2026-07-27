! SPDX-License-Identifier: LGPL-3.0-or-later
program sharper_demo
   use sharper, only: dp, sr_result, sropt_result, test_result
   use sharper, only: fit_sr, fit_sropt, sr_confint, sr_test, sric
   implicit none
   real(dp) :: returns(252,3)
   real(dp), allocatable :: ci(:, :)
   type(sr_result) :: ratios
   type(sropt_result) :: optimum
   type(test_result) :: test
   integer :: i

   do i = 1, size(returns,1)
      returns(i,1) = 0.0005_dp+0.010_dp*sin(0.17_dp*real(i,dp)) + &
                     0.004_dp*cos(0.43_dp*real(i,dp))
      returns(i,2) = 0.0003_dp+0.006_dp*sin(0.11_dp*real(i,dp)) + &
                     0.005_dp*cos(0.37_dp*real(i,dp))
      returns(i,3) = 0.0002_dp+0.004_dp*sin(0.29_dp*real(i,dp)) + &
                     0.007_dp*cos(0.19_dp*real(i,dp))
   end do

   ratios = fit_sr(returns,ope=252.0_dp,higher_order=.true.)
   ci = sr_confint(ratios,0.95_dp,'t')
   optimum = fit_sropt(returns,ope=252.0_dp)
   test = sr_test(fit_sr(returns(:,1),ope=252.0_dp),0.0_dp,'greater','exact')

   print '(a)', 'SharpeR modern Fortran demo'
   print '(a,3f12.5)', 'annualized Sharpe ratios: ',ratios%value
   print '(a,2f12.5)', 'asset 1 95% interval:   ',ci(1,:)
   print '(a,f12.5)', 'optimal Sharpe ratio:    ',optimum%value
   print '(a,f12.5)', 'SRIC:                    ',sric(optimum)
   print '(a,f12.6)', 'asset 1 greater p-value: ',test%p_value
end program sharper_demo
