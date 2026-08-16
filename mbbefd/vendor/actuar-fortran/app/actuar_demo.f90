! SPDX-License-Identifier: GPL-2.0-or-later
program actuar_demo
  use actuar, only : dp, ppareto, qpareto, aggregate_distribution, &
    panjer_poisson, aggregate_var, aggregate_cte
  implicit none
  type(aggregate_distribution) :: aggregate
  real(dp) :: severity(4)

  print '(a,f12.8)', 'Pareto F(10): ',ppareto(10.0_dp,2.5_dp,4.0_dp)
  print '(a,f12.8)', 'Pareto q(0.99): ',qpareto(0.99_dp,2.5_dp,4.0_dp)

  severity=[0.10_dp,0.45_dp,0.30_dp,0.15_dp]
  aggregate=panjer_poisson(severity,2.0_dp,tolerance=0.999999_dp,max_terms=200)
  if(.not.aggregate%ok) error stop trim(aggregate%message)
  print '(a,f10.4)', 'Aggregate mean: ',aggregate%mean()
  print '(a,f10.4)', 'Aggregate VaR 0.95: ',aggregate_var(aggregate,0.95_dp)
  print '(a,f10.4)', 'Aggregate CTE 0.95: ',aggregate_cte(aggregate,0.95_dp)
end program actuar_demo
