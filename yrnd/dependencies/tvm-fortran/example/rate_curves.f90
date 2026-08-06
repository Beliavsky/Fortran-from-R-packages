! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
program rate_curves
   use tvm
   implicit none
   type(rate_curve_t) :: curve

   curve = rate_curve_from_rates([0.10_dp, 0.20_dp, 0.30_dp], "zero_eff")
   print '(a,*(f10.6,1x))', "Discounts: ", curve%discount(curve%knots)
   print '(a,*(f10.6,1x))', "Swap rates: ", curve%rate_grid("swap")
   print '(a,f12.6)', "PV: ", curve%present_value([-1.0_dp, 1.10_dp], [0.0_dp, 1.0_dp])
end program rate_curves
