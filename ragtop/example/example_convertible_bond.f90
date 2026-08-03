program example_convertible_bond
   use ragtop
   implicit none
   type(market_spec) :: market
   type(cashflow_schedule) :: coupons
   type(exercise_schedule) :: calls
   type(instrument_spec) :: bond
   real(dp) :: value

   market%short_rate = 0.04_dp
   market%volatility = 0.35_dp
   market%default_intensity = 0.03_dp
   coupons = make_cashflow_schedule( &
      [0.5_dp,1.0_dp,1.5_dp,2.0_dp,2.5_dp,3.0_dp], &
      [2.0_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp])
   calls = make_exercise_schedule([2.0_dp],[105.0_dp])
   bond = ConvertibleBond(3.0_dp,100.0_dp,0.8_dp,coupons,0.4_dp, &
                          calls=calls,name='Example convertible')
   value = find_present_value(100.0_dp,180,bond,market, &
                              std_devs_width=4.0_dp)
   write(*,'(a,f12.6)') 'Convertible value: ',value
end program example_convertible_bond
