! SPDX-License-Identifier: GPL-2.0-or-later
program vanilla_options
   use fer, only : dp, bachelier_price, black_scholes_price
   implicit none

   real(dp), parameter :: strike = 100.0_dp
   real(dp), parameter :: forward = 105.0_dp
   real(dp), parameter :: texp = 1.0_dp
   real(dp), parameter :: sigma = 0.20_dp
   real(dp), parameter :: df = 0.97_dp

   print '(a,es16.8)', 'Black-Scholes call: ', &
      black_scholes_price(strike, forward, texp, sigma, 1, df)
   print '(a,es16.8)', 'Bachelier call:     ', &
      bachelier_price(strike, forward, texp, 20.0_dp, 1, df)
end program vanilla_options
