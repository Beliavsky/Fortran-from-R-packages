program example_greeks
   use ragtop
   implicit none
   type(market_spec) :: market
   type(instrument_spec) :: option
   type(greek_result) :: greeks

   market%short_rate = 0.03_dp
   market%volatility = 0.25_dp
   market%default_intensity = 0.02_dp
   option = AmericanOption(1.25_dp,100.0_dp,put_option,'Put')
   greeks = find_greeks(95.0_dp,120,option,market)
   write(*,'(a,f12.6)') 'Price:  ',greeks%price
   write(*,'(a,f12.6)') 'Delta:  ',greeks%delta
   write(*,'(a,f12.6)') 'Gamma:  ',greeks%gamma
   write(*,'(a,f12.6)') 'Vega:   ',greeks%vega
   write(*,'(a,f12.6)') 'Rho:    ',greeks%rho
   write(*,'(a,f12.6)') 'Hazard: ',greeks%hazard_sensitivity
end program example_greeks
