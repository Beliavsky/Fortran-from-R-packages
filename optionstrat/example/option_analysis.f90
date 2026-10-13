! SPDX-License-Identifier: GPL-3.0-only
program option_analysis
   use optionstrat
   implicit none
   type(option_result) :: call
   type(vertical_result) :: spread
   type(price_range) :: interval
   real(dp) :: sigma
   call = calleval(100.0_dp, 100.0_dp, 0.2_dp, 0.25_dp, 0.03_dp, 0.01_dp)
   sigma = implied_volatility('call', call%premium, 100.0_dp, 100.0_dp, 0.25_dp, 0.03_dp, 0.01_dp)
   print '(a,f12.6)', 'European call premium: ', call%premium
   print '(a,f12.6)', 'Delta:                 ', call%delta
   print '(a,f12.6)', 'Recovered volatility:  ', sigma
   print '(a,f12.6)', 'R-compatible theta:    ', call%theta
   print '(a,f12.6)', 'Standard daily theta:  ', &
      calltheta_standard(100.0_dp, 100.0_dp, 0.2_dp, 0.25_dp, 0.03_dp, 0.01_dp)
   spread = vertical('call', 100.0_dp, 110.0_dp, 90.0_dp, 0.25_dp, 0.03_dp, 0.2_dp)
   print '(a,f12.6)', 'Bull-call maximum P/L: ', spread%max_profit
   print '(a,f12.6)', 'Bull-call minimum P/L: ', spread%max_loss
   interval = prob_btwn_quantile_standard(100.0_dp, 0.8_dp, asd = 0.2_dp, dte = 30.0_dp, mean = 0.01_dp)
   print '(a,2f12.6)', '80% arithmetic-normal price interval: ', interval%lower_price, interval%upper_price
end program option_analysis
