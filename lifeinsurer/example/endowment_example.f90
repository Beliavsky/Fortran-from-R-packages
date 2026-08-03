! SPDX-License-Identifier: GPL-2.0-or-later
program endowment_example
  use lifeinsurer
  implicit none
  type(insurance_tariff)::t
  type(mortality_table)::mt
  type(contract_result)::r
  integer::i
  allocate(mt%qx(20)); mt%qx=[(0.002_dp+0.0004_dp*real(i-1,dp),i=1,20)]
  t%tariff_type=tariff_endowment; t%policy_period=20; t%premium_period=20
  t%sum_insured=100000.0_dp; t%interest=0.02_dp
  t%costs=initialize_costs(alpha=0.035_dp,beta=0.02_dp,gamma_contract=0.0005_dp,unit_cost=12.0_dp,security=0.03_dp)
  call calculate_contract(t,mt,r)
  print '(a,f12.2)','Written annual premium: ',r%premiums%written_yearly
  print '(a,f12.2)','Contractual reserve after 10 years: ',r%reserves%contractual(11)
end program
