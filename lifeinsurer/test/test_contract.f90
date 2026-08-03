! SPDX-License-Identifier: GPL-2.0-or-later
program test_contract
  use lifeinsurer
  implicit none
  type(insurance_tariff)::t,u
  type(mortality_table)::mt
  type(contract_result)::r
  type(profit_rate_table)::rates
  integer::n
  n=20; allocate(mt%qx(n)); mt%qx=0.01_dp
  t%tariff_type=tariff_endowment; t%policy_period=n; t%premium_period=n
  t%sum_insured=10000.0_dp; t%interest=0.03_dp
  t%costs=initialize_costs(alpha=0.04_dp,zillmer=0.02_dp,beta=0.01_dp,gamma_contract=0.0005_dp,unit_cost=10.0_dp,security=0.02_dp)
  t%rounding%premium_gross=2; t%rounding%premium_net=2; t%rounding%premium_written=2
  call calculate_contract(t,mt,r)
  call check(r%status==0,'contract status')
  call check(r%premiums%gross>r%premiums%net,'gross exceeds net')
  call check(abs(r%reserves%net(n+1)-t%sum_insured*(1.0_dp+t%costs%security))<1e-6_dp, &
    'maturity reserve before payment')
  call check(r%reserves%contractual(1)>=0.0_dp,'nonnegative reserve')
  u=premium_waiver(t); call check(u%premium_period==0 .and. u%premium_waiver,'waiver')
  u=extend_contract(t,5,new_sum_insured=15000.0_dp,new_interest=0.01_dp)
  call check(u%policy_period==25 .and. abs(u%interest-0.01_dp)<1e-14_dp,'extension')
  allocate(rates%guaranteed_interest(n+1),rates%interest_profit(n+1),rates%interest_profit2(n+1))
  allocate(rates%mortality_profit(n+1),rates%expense_profit(n+1),rates%sum_profit(n+1))
  allocate(rates%terminal_bonus(n+1),rates%terminal_bonus_fund(n+1))
  rates%guaranteed_interest=0.01_dp; rates%interest_profit=0.005_dp; rates%interest_profit2=0.0_dp
  rates%mortality_profit=0.02_dp; rates%expense_profit=0.001_dp; rates%sum_profit=0.0_dp
  rates%terminal_bonus=0.1_dp; rates%terminal_bonus_fund=0.05_dp
  call calculate_contract(t,mt,r,rates)
  call check(r%status==0 .and. all(r%profits%total_assignment>=0.0_dp),'profit participation')
  print '(a)','test_contract: PASS'
contains
  subroutine check(ok,msg); logical,intent(in)::ok; character(*),intent(in)::msg
    if(.not.ok) then; print '(a,1x,a)','FAIL:',msg; error stop 1; end if
  end subroutine
end program
