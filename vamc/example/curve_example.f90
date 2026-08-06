program curve_example
  use vamc
  implicit none
  real(dp) :: rates(8)
  integer :: tenors(8), i
  type(date_type) :: curve_date
  type(yield_curve_type) :: curve
  type(status_type) :: status
  rates = [0.69_dp,0.77_dp,0.88_dp,1.01_dp,1.14_dp,1.38_dp,1.66_dp,2.15_dp]*0.01_dp
  tenors = [1,2,3,4,5,7,10,30]
  call parse_date('2016-02-08',curve_date,status)
  call build_curve_named(rates,tenors,6,'Thirty360',6,'ACT360','NY','Modified_Foll',curve_date,2, &
                         'Thirty360',curve,status=status)
  if (.not. status%ok()) error stop status%message
  do i=1,size(curve%forward_rates)
    write(*,'(a,1x,f12.9)') curve%observation_dates(i)%to_string(),curve%forward_rates(i)
  end do
end program curve_example
