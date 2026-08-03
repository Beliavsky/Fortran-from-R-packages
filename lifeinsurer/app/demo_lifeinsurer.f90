! SPDX-License-Identifier: GPL-2.0-or-later
program demo_lifeinsurer
  use lifeinsurer
  implicit none
  type(insurance_tariff)::t
  type(mortality_table)::mt
  type(contract_result)::r
  integer::i
  allocate(mt%qx(30)); mt%qx=[(0.0015_dp*1.08_dp**real(i-1,dp),i=1,30)]
  t%tariff_type=tariff_endowment; t%policy_period=30; t%premium_period=25
  t%sum_insured=100000.0_dp; t%interest=0.02_dp; t%premium_refund=1.0_dp
  t%costs=initialize_costs(alpha=0.04_dp,zillmer=0.025_dp,beta=0.015_dp,gamma_contract=0.0006_dp,unit_cost=15.0_dp,security=0.03_dp)
  call calculate_contract(t,mt,r)
  print '(a,i0)','Status: ',r%status
  print '(a,f12.2)','Net premium: ',r%premiums%net
  print '(a,f12.2)','Gross premium: ',r%premiums%gross
  print '(a,f12.2)','Written installment: ',r%premiums%written
  print '(a,f12.2)','Contractual reserve at year 15: ',r%reserves%contractual(16)
end program
