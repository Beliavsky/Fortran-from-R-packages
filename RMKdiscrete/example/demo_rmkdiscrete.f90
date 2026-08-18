! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of RMKdiscrete 0.1 by Robert M. Kirkpatrick.
! See LICENSE, COPYING, and TRANSLATION_NOTES.md for provenance.
program demo_rmkdiscrete
  use rmkdiscrete, only : dp, dlgp, plgp, qlgp, lgp_summary, slgp, dbilgp
  implicit none
  type(lgp_summary) :: s
  real(dp) :: theta(3),lambda(3)

  print '(a,f12.8)','P(X=3), LGP(theta=2, lambda=0.3): ',dlgp(3,2.0_dp,0.3_dp)
  print '(a,f12.8)','P(X<=3): ',plgp(3.0_dp,2.0_dp,0.3_dp)
  print '(a,f8.1)','median: ',qlgp(0.5_dp,2.0_dp,0.3_dp)
  s=slgp(2.0_dp,0.3_dp)
  print '(a,f12.6,a,f12.6)','mean=',s%mean,' variance=',s%variance

  theta=[1.0_dp,2.0_dp,1.5_dp]
  lambda=[0.2_dp,0.1_dp,0.0_dp]
  print '(a,f12.8)','P(Y1=2,Y2=1), bivariate LGP: ',dbilgp(2,1,theta,lambda)
end program demo_rmkdiscrete
