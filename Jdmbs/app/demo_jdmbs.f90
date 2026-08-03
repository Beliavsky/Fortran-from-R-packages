! SPDX-License-Identifier: GPL-2.0-or-later
program demo_jdmbs
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: diffusion, independent_jump, network_jump
   real(dp) :: transmission(2, 2)
   real(dp) :: spot(2), mu(2), sigma(2), strike(2)

   spot = [100.0_dp, 120.0_dp]
   mu = [0.04_dp, 0.03_dp]
   sigma = [0.20_dp, 0.24_dp]
   strike = [100.0_dp, 125.0_dp]
   transmission = reshape([1.0_dp, 0.5_dp, 0.4_dp, 1.0_dp], [2, 2])
   control%day = 90
   control%monte_carlo = 10000
   control%seed = 123456789_int64
   control%discount_rate = 0.03_dp
   control%legacy_mode = .false.

   call normal_bs(spot, mu, sigma, strike, diffusion, control)
   call jdm_bs(spot, mu, sigma, 2.0_dp, strike, independent_jump, control)
   call jdm_new_bs(transmission, spot, mu, sigma, 2.0_dp, strike, &
      network_jump, control)

   write (*, '(a)') 'Jdmbs modern Fortran demo'
   write (*, '(a,2f11.4)') 'diffusion calls:       ', diffusion%call_price
   write (*, '(a,2f11.4)') 'independent-jump calls:', independent_jump%call_price
   write (*, '(a,2f11.4)') 'network-jump calls:    ', network_jump%call_price
end program demo_jdmbs
