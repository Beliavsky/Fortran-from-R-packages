program test_portfolio_aging
  use vamc
  use vamc_test_support
  implicit none
  type(date_type) :: birth_range(2), issue_range(2), dates(4), target
  type(portfolio_type) :: portfolio, aged
  type(mortality_table_type) :: table
  type(status_type) :: status
  character(len=4) :: products(2)
  real(dp) :: percentages(2), rider(2), rollup(2), withdrawal(2), fund_fees(2)
  real(dp) :: historical(4,2), discount(24)
  integer :: i
  birth_range=[make_date(1960,1,1),make_date(1980,12,31)]
  issue_range=[make_date(2020,1,1),make_date(2020,1,31)]
  products=['DBRP','WBSU']
  percentages=[0.5_dp,0.5_dp]
  rider=[25.0_dp,50.0_dp]
  rollup=5.0_dp
  withdrawal=5.0_dp
  fund_fees=[30.0_dp,40.0_dp]
  call gen_port_inception(birth_range,issue_range,[1,1],[100.0_dp,100.0_dp],0.5_dp,fund_fees,200.0_dp, &
                          percentages,products,rider,rollup,withdrawal,1,portfolio,seed=7,status=status)
  call assert_true(status%ok(),'portfolio generation')
  call assert_true(portfolio%size()==2,'one policy per product type')
  do i=1,portfolio%size()
    call assert_close(portfolio%policies(i)%account_value(),100.0_dp,1.0e-13_dp,'generated account value')
  end do
  call make_test_mortality(table)
  do i=1,4
    dates(i)=make_date(2020,i,1)
  end do
  target=make_date(2020,4,1)
  historical=1.0_dp
  discount=1.0_dp
  call age_portfolio(portfolio,table,historical,dates,1.0_dp/12.0_dp,target,discount,aged,status=status)
  call assert_true(status%ok(),'portfolio aging')
  call assert_true(aged%policies(1)%current_date>portfolio%policies(1)%current_date,'current date advanced')
  call assert_true(aged%policies(1)%account_value()<portfolio%policies(1)%account_value(),'fees reduce account value')
  print '(a)','test_portfolio_aging: PASS'
end program test_portfolio_aging
