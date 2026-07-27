! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program bcc1997_demo
   use bcc1997
   implicit none

   type(bcc_parameters) :: parameters
   type(bcc_result) :: result

   parameters = bcc_parameters(kappa_v=1.5_dp, kappa_r=0.4_dp, &
      theta_v=0.04_dp, theta_r=0.03_dp, sigma_v=0.3_dp, sigma_r=0.1_dp, &
      mu_j=-0.05_dp, sigma_j=0.2_dp, rho=-0.6_dp, lambda=0.2_dp, &
      spot=100.0_dp, strike=105.0_dp, variance0=0.04_dp, &
      rate0=0.025_dp, maturity=1.25_dp)

   result = bcc_price(parameters)
   if (result%status /= 0) then
      print '(a)', trim(result%message)
      error stop 1
   end if

   print '(a,f12.6)', 'call price: ', result%call
   print '(a,f12.6)', 'put price:  ', result%put
   print '(a,f12.8)', 'probability 1: ', result%probability1
   print '(a,f12.8)', 'probability 2: ', result%probability2
   print '(a,i0)', 'integrand evaluations: ', &
      result%evaluations1 + result%evaluations2
end program bcc1997_demo
