! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program vares_demo
  use vares, only : dp, varnormal, esnormal, varpareto, espareto, varweibull, esweibull
  implicit none
  real(dp), parameter :: p = 0.05_dp
  print '(a,f14.8)', 'Normal lower-tail VaR: ', varnormal(p)
  print '(a,f14.8)', 'Normal lower-tail ES:  ', esnormal(p)
  print '(a,f14.8)', 'Pareto VaR:            ', varpareto(p, 1.0_dp, 3.0_dp)
  print '(a,f14.8)', 'Pareto ES:             ', espareto(p, 1.0_dp, 3.0_dp)
  print '(a,f14.8)', 'Weibull VaR:           ', varweibull(p, 2.0_dp, 1.5_dp)
  print '(a,f14.8)', 'Weibull ES:            ', esweibull(p, 2.0_dp, 1.5_dp)
end program vares_demo
