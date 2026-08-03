! SPDX-License-Identifier: GPL-2.0-or-later
program jump_diffusion_example
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: price
   real(dp) :: spot(3), drift(3), volatility(3), strike(3)

   spot = [100.0_dp, 120.0_dp, 75.0_dp]
   drift = [0.03_dp, 0.02_dp, 0.04_dp]
   volatility = [0.18_dp, 0.22_dp, 0.28_dp]
   strike = [100.0_dp, 125.0_dp, 80.0_dp]
   control%day = 180
   control%monte_carlo = 15000
   control%seed = 314159265_int64
   control%discount_rate = 0.03_dp
   control%legacy_mode = .false.

   call jdm_bs(spot, drift, volatility, 2.0_dp, strike, price, control)
   write (*, '(a,3f12.5)') 'jump call prices: ', price%call_price
   write (*, '(a,3f12.5)') 'jump put prices:  ', price%put_price
   write (*, '(a,i0)') 'simulated jump events: ', price%jump_events
end program jump_diffusion_example
