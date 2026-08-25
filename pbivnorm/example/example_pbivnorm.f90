! SPDX-License-Identifier: GPL-2.0-or-later
program example_pbivnorm
   use pbivnorm_mod, only : dp, pbivnorm
   implicit none
   real(dp) :: x(4), y(4), p(4)

   x = [-1.0_dp, 0.0_dp, 0.5_dp, 1.0_dp]
   y = [ 0.5_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   p = pbivnorm(x, y, 0.5_dp)

   print '(a)', '      x       y      rho       CDF'
   print '(3f8.3,f12.8)', x(1), y(1), 0.5_dp, p(1)
   print '(3f8.3,f12.8)', x(2), y(2), 0.5_dp, p(2)
   print '(3f8.3,f12.8)', x(3), y(3), 0.5_dp, p(3)
   print '(3f8.3,f12.8)', x(4), y(4), 0.5_dp, p(4)
end program example_pbivnorm
