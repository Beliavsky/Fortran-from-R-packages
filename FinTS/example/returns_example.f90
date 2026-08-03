! SPDX-License-Identifier: GPL-2.0-or-later
program returns_example
   use fints
   implicit none
   real(dp) :: r(3)
   r = [0.01_dp, -0.02_dp, 0.03_dp]
   print '(3f12.8)', simple2logReturns(r)
end program returns_example
