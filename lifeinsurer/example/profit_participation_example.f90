! SPDX-License-Identifier: GPL-2.0-or-later
program profit_participation_example
  use lifeinsurer
  implicit none
  type(insurance_tariff)::t
  type(mortality_table)::mt
  type(contract_result)::r
  type(profit_rate_table)::pr
  integer::n
  n=15; allocate(mt%qx(n)); mt%qx=0.005_dp
  t%tariff_type=tariff_endowment; t%policy_period=n; t%premium_period=n; t%sum_insured=50000.0_dp; t%interest=0.015_dp
  allocate(pr%guaranteed_interest(n+1),pr%interest_profit(n+1),pr%interest_profit2(n+1))
  allocate(pr%mortality_profit(n+1),pr%expense_profit(n+1),pr%sum_profit(n+1))
  allocate(pr%terminal_bonus(n+1),pr%terminal_bonus_fund(n+1))
  pr%guaranteed_interest=0.015_dp; pr%interest_profit=0.006_dp; pr%interest_profit2=0.0_dp
  pr%mortality_profit=0.05_dp; pr%expense_profit=0.001_dp; pr%sum_profit=0.0_dp
  pr%terminal_bonus=0.1_dp; pr%terminal_bonus_fund=0.05_dp
  call calculate_contract(t,mt,r,pr)
  print '(a,f12.2)','Profit benefit at maturity: ',r%profits%benefit(n+1)
end program
