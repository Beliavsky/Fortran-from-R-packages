! SPDX-License-Identifier: GPL-2.0-or-later
program normal_option_example
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: price
   real(dp) :: spot(2), drift(2), volatility(2), strike(2)

   spot = [100.0_dp, 80.0_dp]
   drift = [0.04_dp, 0.02_dp]
   volatility = [0.20_dp, 0.30_dp]
   strike = [100.0_dp, 85.0_dp]
   control%day = 90
   control%monte_carlo = 20000
   control%seed = 20260731_int64
   control%discount_rate = 0.04_dp
   control%legacy_mode = .false.

   call normal_bs(spot, drift, volatility, strike, price, control)
   write (*, '(a,2f12.5)') 'call prices: ', price%call_price
   write (*, '(a,2f12.5)') 'put prices:  ', price%put_price
   write (*, '(a,2f12.5)') 'call SE:     ', price%call_se
end program normal_option_example
