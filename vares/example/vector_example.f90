! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
program vector_example
  use vares, only : dp, varnormal, pnormal
  implicit none
  real(dp) :: p(5), q(5), back(5)
  p = [0.01_dp, 0.05_dp, 0.25_dp, 0.50_dp, 0.95_dp]
  q = varnormal(p, mu=1.0_dp, sigma=2.0_dp)
  back = pnormal(q, mu=1.0_dp, sigma=2.0_dp)
  print '(3f16.8)', p, q, back
end program vector_example
