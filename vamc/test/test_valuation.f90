program test_valuation
  use vamc
  use vamc_test_support
  implicit none
  type(policy_type) :: policy
  type(portfolio_type) :: portfolio
  type(mortality_table_type) :: table
  type(policy_valuation_type) :: v2, v3
  type(portfolio_valuation_type) :: vp
  type(status_type) :: status
  real(dp) :: scenario2(12,2), scenario3(2,12,2), discount(12)
  call make_test_policy('DBRP',policy)
  call make_test_mortality(table)
  scenario2=1.0_dp
  scenario3(1,:,:)=scenario2
  scenario3(2,:,:)=scenario2
  discount=1.0_dp
  call valuate_one_policy(policy,table,scenario2,1.0_dp/12.0_dp,discount,v2,status=status)
  call assert_true(status%ok(),'single-scenario valuation')
  call valuate_one_policy(policy,table,scenario3,1.0_dp/12.0_dp,discount,v3,status=status)
  call assert_true(status%ok(),'multi-scenario valuation')
  call assert_close(v2%policy_value,v3%policy_value,1.0e-13_dp,'scenario averaging')
  allocate(portfolio%policies(2))
  portfolio%policies(1)=policy
  portfolio%policies(2)=policy
  portfolio%policies(2)%record_id=2
  call valuate_portfolio(portfolio,table,scenario3,1.0_dp/12.0_dp,discount,vp,status=status)
  call assert_true(status%ok(),'portfolio valuation')
  call assert_close(vp%portfolio_value,2.0_dp*v3%policy_value,1.0e-13_dp,'portfolio aggregation')
  call assert_close(vp%portfolio_risk_charge,2.0_dp*v3%risk_charge,1.0e-13_dp,'risk-charge aggregation')
  print '(a)','test_valuation: PASS'
end program test_valuation
