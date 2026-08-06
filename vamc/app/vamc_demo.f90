program vamc_demo
  use vamc
  implicit none
  real(dp) :: rates(8), covariance(2,2), forward(24), fund_map(2,2), discount(24)
  real(dp), allocatable :: index_scenarios(:,:,:), fund_scenarios(:,:,:)
  integer :: tenors(8), ages(121), i
  real(dp) :: female_q(121), male_q(121)
  type(date_type) :: curve_date
  type(yield_curve_type) :: curve
  type(policy_type) :: policy
  type(mortality_table_type) :: mortality
  type(policy_valuation_type) :: value
  type(status_type) :: status

  rates = [0.69_dp,0.77_dp,0.88_dp,1.01_dp,1.14_dp,1.38_dp,1.66_dp,2.15_dp]*0.01_dp
  tenors = [1,2,3,4,5,7,10,30]
  call parse_date('2016-02-08',curve_date,status)
  call build_curve_named(rates,tenors,6,'Thirty360',6,'ACT360','NY','Modified_Foll',curve_date,2, &
                         'Thirty360',curve,status=status)
  if (.not. status%ok()) error stop status%message

  covariance(1,:) = [0.0225_dp,0.0060_dp]
  covariance(2,:) = [0.0060_dp,0.0400_dp]
  forward = 0.02_dp
  call gen_index_scen(covariance,100,24,1.0_dp/12.0_dp,forward,index_scenarios,seed=42,status=status)
  if (.not. status%ok()) error stop status%message
  fund_map(1,:) = [1.0_dp,0.0_dp]
  fund_map(2,:) = [0.4_dp,0.6_dp]
  call gen_fund_scen(fund_map,index_scenarios,fund_scenarios,status)
  if (.not. status%ok()) error stop status%message

  policy%record_id = 1
  policy%gender = 'F'
  policy%product_type = 'DBRP'
  policy%birth_date = make_date(1970,1,1)
  policy%issue_date = make_date(2020,1,1)
  policy%current_date = make_date(2020,1,1)
  policy%maturity_date = make_date(2022,1,1)
  policy%rider_fee = 0.0025_dp
  policy%base_fee = 0.0200_dp
  policy%guarantee_amount = 120000.0_dp
  allocate(policy%fund_numbers(2),policy%fund_values(2),policy%fund_fees(2))
  policy%fund_numbers = [1,2]
  policy%fund_values = [60000.0_dp,40000.0_dp]
  policy%fund_fees = [0.0030_dp,0.0040_dp]

  do i=1,121
    ages(i)=i-1
    female_q(i)=min(1.0_dp,0.0005_dp*exp(0.075_dp*real(max(0,ages(i)-30),dp)))
    male_q(i)=min(1.0_dp,1.2_dp*female_q(i))
  end do
  female_q(121)=1.0_dp
  male_q(121)=1.0_dp
  call make_mortality_table(ages,female_q,male_q,mortality,status)
  discount = [(exp(-0.02_dp*real(i,dp)/12.0_dp),i=1,24)]
  call valuate_one_policy(policy,mortality,fund_scenarios,1.0_dp/12.0_dp,discount,value,status=status)
  if (.not. status%ok()) error stop status%message

  write(*,'(a,9f12.8)') 'Forward curve: ',curve%forward_rates
  write(*,'(a,f14.2)') 'Estimated policy value: ',value%policy_value
  write(*,'(a,f14.2)') 'Estimated risk charge: ',value%risk_charge
end program vamc_demo
