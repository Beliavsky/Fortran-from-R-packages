! SPDX-License-Identifier: GPL-2.0-or-later
program correlated_companies_example
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: price
   real(dp) :: transmission(3, 3)
   real(dp) :: spot(3), drift(3), volatility(3), strike(3)

   transmission = reshape([ &
      1.00_dp, 0.65_dp, 0.25_dp, &
      0.55_dp, 1.00_dp, 0.40_dp, &
      0.20_dp, 0.35_dp, 1.00_dp], [3, 3])
   spot = [100.0_dp, 95.0_dp, 110.0_dp]
   drift = 0.03_dp
   volatility = [0.20_dp, 0.18_dp, 0.25_dp]
   strike = [100.0_dp, 100.0_dp, 105.0_dp]
   control%day = 120
   control%monte_carlo = 12000
   control%seed = 271828182_int64
   control%discount_rate = 0.03_dp
   control%legacy_mode = .false.

   call jdm_new_bs(transmission, spot, drift, volatility, 3.0_dp, &
      strike, price, control)
   write (*, '(a,3f12.5)') 'network-jump calls: ', price%call_price
   write (*, '(a,3f12.5)') 'network-jump puts:  ', price%put_price
end program correlated_companies_example
