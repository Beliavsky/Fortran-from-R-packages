! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program strike_curve
   use bcc1997
   implicit none

   type(bcc_parameters) :: parameters
   type(bcc_result), allocatable :: results(:)
   real(dp), parameter :: strikes(7) = [70.0_dp, 80.0_dp, 90.0_dp, &
      100.0_dp, 110.0_dp, 120.0_dp, 130.0_dp]
   integer :: i

   parameters = bcc_parameters(kappa_v=1.4_dp, kappa_r=0.3_dp, &
      theta_v=0.04_dp, theta_r=0.025_dp, sigma_v=0.28_dp, &
      sigma_r=0.08_dp, mu_j=-0.04_dp, sigma_j=0.18_dp, rho=-0.55_dp, &
      lambda=0.18_dp, spot=100.0_dp, strike=100.0_dp, &
      variance0=0.04_dp, rate0=0.02_dp, maturity=1.0_dp)

   call bcc_price_strikes(parameters, strikes, results)
   print '(a)', ' strike         call          put'
   do i = 1, size(strikes)
      print '(f7.2,2f13.6)', strikes(i), results(i)%call, results(i)%put
   end do
end program strike_curve
