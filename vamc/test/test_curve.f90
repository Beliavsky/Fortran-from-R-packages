program test_curve
  use vamc
  use vamc_test_support
  implicit none
  real(dp) :: rates(8), expected(9)
  integer :: tenors(8)
  type(date_type) :: date, rolled
  type(yield_curve_type) :: curve
  type(status_type) :: status
  rates = [0.69_dp,0.77_dp,0.88_dp,1.01_dp,1.14_dp,1.38_dp,1.66_dp,2.15_dp]*0.01_dp
  tenors = [1,2,3,4,5,7,10,30]
  expected = [0.006912035_dp,0.008520036_dp,0.011060713_dp,0.014146403_dp,0.016846310_dp, &
              0.020451281_dp,0.024514485_dp,0.032098279_dp,0.032098279_dp]
  call parse_date('2016-02-08',date,status)
  call assert_true(status%ok(),'parse curve date')
  call build_curve_named(rates,tenors,6,'Thirty360',6,'ACT360','NY','Modified_Foll',date,2,'Thirty360',curve,status=status)
  call assert_true(status%ok(),'build curve')
  call assert_all_close(curve%forward_rates,expected,1.0e-7_dp,'published forward curve fixture')
  call assert_close(frac_year(make_date(2020,1,1),make_date(2021,1,1),dcc_actact),1.0_dp,1.0e-14_dp,'ACT/ACT leap year')
  rolled=roll_date(make_date(2021,7,4),bdc_following,calendar_ny)
  call assert_true(rolled==make_date(2021,7,6),'NY following holiday roll')
  print '(a)','test_curve: PASS'
end program test_curve
