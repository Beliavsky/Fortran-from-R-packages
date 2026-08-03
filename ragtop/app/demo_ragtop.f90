program demo_ragtop
   use ragtop
   implicit none
   type(market_spec) :: market
   type(cashflow_schedule) :: coupons
   type(instrument_spec) :: euro, amer, convertible
   type(option_value) :: bs
   type(greek_result) :: greeks
   real(dp) :: pde_euro, amer_price, convertible_price, implied

   market%short_rate = 0.035_dp
   market%volatility = 0.30_dp
   market%default_intensity = 0.025_dp
   market%dividend_rate = 0.01_dp

   euro = EuropeanOption(1.5_dp,100.0_dp,call_option,'European call')
   amer = AmericanOption(1.5_dp,100.0_dp,put_option,'American put')
   coupons = make_cashflow_schedule( &
      [0.5_dp,1.0_dp,1.5_dp,2.0_dp,2.5_dp,3.0_dp], &
      [2.0_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp])
   convertible = ConvertibleBond(3.0_dp,100.0_dp,0.85_dp,coupons,0.40_dp, &
                                 name='Convertible')

   bs = black_scholes_on_term_structures(call_option,100.0_dp,100.0_dp, &
                                          1.5_dp,market)
   pde_euro = find_present_value(100.0_dp,140,euro,market)
   amer_price = american(put_option,100.0_dp,100.0_dp,1.5_dp,market,140)
   convertible_price = find_present_value(100.0_dp,180,convertible,market)
   implied = implied_volatility_with_term_struct(bs%price,call_option, &
                                                 100.0_dp,100.0_dp,1.5_dp, &
                                                 market)
   greeks = find_greeks(100.0_dp,100,amer,market)

   write(*,'(a)') 'ragtop-fortran demonstration'
   write(*,'(a,f12.6)') 'Black-Scholes European call: ',bs%price
   write(*,'(a,f12.6)') 'PDE European call:           ',pde_euro
   write(*,'(a,f12.6)') 'American put:                ',amer_price
   write(*,'(a,f12.6)') 'Convertible bond:            ',convertible_price
   write(*,'(a,f12.6)') 'Recovered implied vol:       ',implied
   write(*,'(a,f12.6)') 'American put delta:          ',greeks%delta
end program demo_ragtop
