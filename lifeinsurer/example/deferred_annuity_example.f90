! SPDX-License-Identifier: GPL-2.0-or-later
program deferred_annuity_example
  use lifeinsurer
  implicit none
  type(insurance_tariff)::t
  type(mortality_table)::mt
  type(contract_result)::r
  integer::i
  allocate(mt%qx(80)); mt%qx=[(min(0.35_dp,0.001_dp*1.07_dp**real(i-1,dp)),i=1,80)]
  t%tariff_type=tariff_annuity; t%policy_period=80; t%deferral_period=25
  t%premium_period=25; t%premium_refund=1.0_dp; t%sum_insured=1200.0_dp; t%interest=0.005_dp
  call calculate_contract(t,mt,r)
  print '(a,f12.2)','Annual annuity unit gross premium: ',r%premiums%gross
  print '(a,f12.6)','PV of annuity benefits: ',r%present_values%benefits(1)
end program
