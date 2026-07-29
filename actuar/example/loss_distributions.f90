! SPDX-License-Identifier: GPL-2.0-or-later
program loss_distributions
  use actuar, only : dp, dburr, pburr, qburr, mburr, &
    dinvgamma, pinvgamma, qinvgamma
  implicit none
  print '(a,es14.6)', 'Burr density: ',dburr(5.0_dp,2.0_dp,1.5_dp,3.0_dp)
  print '(a,f12.8)', 'Burr CDF: ',pburr(5.0_dp,2.0_dp,1.5_dp,3.0_dp)
  print '(a,f12.6)', 'Burr q(0.95): ',qburr(0.95_dp,2.0_dp,1.5_dp,3.0_dp)
  print '(a,f12.6)', 'Burr mean: ',mburr(1.0_dp,2.0_dp,1.5_dp,3.0_dp)
  print '(a,es14.6)', 'Inverse gamma density: ',dinvgamma(2.0_dp,4.0_dp,3.0_dp)
  print '(a,f12.8)', 'Inverse gamma CDF: ',pinvgamma(2.0_dp,4.0_dp,3.0_dp)
  print '(a,f12.6)', 'Inverse gamma q(0.95): ',qinvgamma(0.95_dp,4.0_dp,3.0_dp)
end program loss_distributions
