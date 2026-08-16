! SPDX-License-Identifier: GPL-2.0-or-later
program aggregate_and_credibility
  use actuar, only : dp, aggregate_distribution, credibility_result, &
    panjer_poisson, aggregate_var, buhlmann_straub
  implicit none
  type(aggregate_distribution) :: aggregate
  type(credibility_result) :: cred
  real(dp) :: severity(3),ratios(3,4),weights(3,4)

  severity=[0.25_dp,0.50_dp,0.25_dp]
  aggregate=panjer_poisson(severity,1.8_dp,tolerance=0.999999_dp,max_terms=150)
  print '(a,f10.4)', 'Aggregate 99% VaR: ',aggregate_var(aggregate,0.99_dp)

  ratios=reshape([1.0_dp,1.1_dp,0.9_dp,1.2_dp,1.3_dp,1.1_dp, &
    0.8_dp,0.9_dp,1.0_dp,1.1_dp,1.2_dp,1.0_dp],[3,4])
  weights=1.0_dp
  cred=buhlmann_straub(ratios,weights)
  if(.not.cred%ok) error stop trim(cred%message)
  print '(a,f10.6)', 'Collective mean: ',cred%collective_mean
  print '(a,3f10.6)', 'Credibility premiums: ',cred%estimates
end program aggregate_and_credibility
