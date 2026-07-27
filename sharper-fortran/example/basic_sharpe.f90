! SPDX-License-Identifier: LGPL-3.0-or-later
program basic_sharpe
   use sharper, only: dp, sr_result, test_result
   use sharper, only: fit_sr, sr_standard_error, sr_test
   implicit none
   real(dp), parameter :: returns(12) = [0.01_dp,0.02_dp,-0.01_dp,0.03_dp, &
      0.0_dp,0.015_dp,0.005_dp,0.025_dp,-0.004_dp,0.018_dp,0.011_dp,0.007_dp]
   real(dp), allocatable :: se(:)
   type(sr_result) :: z
   type(test_result) :: test

   z = fit_sr(returns,ope=12.0_dp,higher_order=.true.)
   se = sr_standard_error(z,'mertens')
   test = sr_test(z,0.0_dp,'greater','exact')
   print '(a,f10.5)', 'annualized Sharpe ratio: ',z%value(1)
   print '(a,f10.5)', 'Mertens standard error:  ',se(1)
   print '(a,f10.6)', 'one-sided p-value:       ',test%p_value
end program basic_sharpe
