program example_term_structures
   use ragtop
   implicit none
   type(market_spec) :: market
   type(option_value) :: value
   integer :: status

   call initialize_discount_curve(market%rates, &
      [0.5_dp,1.0_dp,2.0_dp,5.0_dp], &
      [0.025_dp,0.027_dp,0.030_dp,0.035_dp],status)
   call initialize_volatility_curve(market%vols, &
      [0.5_dp,1.0_dp,2.0_dp,5.0_dp], &
      [0.22_dp,0.24_dp,0.27_dp,0.30_dp],status)
   market%use_rate_curve = .true.
   market%use_vol_curve = .true.
   market%default_intensity = 0.015_dp

   value = black_scholes_on_term_structures(call_option,100.0_dp,100.0_dp, &
                                             2.0_dp,market)
   write(*,'(a,f12.6)') 'Term-structure call value: ',value%price
end program example_term_structures
