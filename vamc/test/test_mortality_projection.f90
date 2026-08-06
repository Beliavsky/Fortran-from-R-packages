program test_mortality_projection
  use vamc
  use vamc_test_support
  implicit none
  type(policy_type) :: policy
  type(mortality_table_type) :: table
  type(mortality_factors_type) :: mortality
  type(projection_result_type) :: result
  type(valuation_options_type) :: corrected
  type(status_type) :: status
  real(dp) :: scenario(12,2), pq(12), p(12), df(12)
  character(len=4) :: products(19)
  integer :: i
  products=default_product_types()
  call make_test_mortality(table)
  call make_test_policy('DBRP',policy)
  call calc_mort_factors(policy,table,1.0_dp/12.0_dp,mortality,status)
  call assert_true(status%ok(),'mortality factors')
  call assert_true(size(mortality%pq)>=12 .and. all(mortality%pq>=0.0_dp),'mortality horizon and positivity')
  scenario=1.0_dp
  pq=1.0_dp
  p=1.0_dp
  df=1.0_dp
  call project_policy(policy,scenario(1:2,:),1.0_dp/12.0_dp,pq,p,df,result,status=status)
  call assert_true(status%ok(),'DBRP projection')
  call assert_close(result%death_benefit,40.0_dp,1.0e-12_dp,'DBRP death benefit')
  call make_test_policy('DBRU',policy)
  call project_policy(policy,scenario(1:2,:),1.0_dp/12.0_dp,pq,p,df,result,status=status)
  call assert_close(result%death_benefit,52.0_dp,1.0e-12_dp,'DBRU annual roll-up')
  corrected=valuation_options_type(source_compatible_timing=.false.,source_compatible_ab_renewal=.false., &
                                    source_compatible_maturity_vector=.false.,source_compatible_aging_index=.false.)
  call make_test_policy('MBRP',policy)
  call project_policy(policy,scenario,1.0_dp/12.0_dp,pq,p,df,result,corrected,status)
  call assert_close(result%living_benefit,20.0_dp,1.0e-12_dp,'corrected maturity benefit')
  do i=1,size(products)
    call make_test_policy(products(i),policy)
    call project_policy(policy,scenario,1.0_dp/12.0_dp,pq,p,df,result,status=status)
    call assert_true(status%ok(),'all rider types dispatch')
    call assert_true(result%death_benefit>=0.0_dp .and. result%living_benefit>=0.0_dp .and. &
                     result%risk_charge>=0.0_dp,'nonnegative rider outputs')
  end do
  print '(a)','test_mortality_projection: PASS'
end program test_mortality_projection
