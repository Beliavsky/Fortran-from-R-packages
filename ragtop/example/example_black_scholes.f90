program example_black_scholes
   use ragtop
   implicit none
   type(option_value) :: value
   real(dp) :: implied

   value = black_scholes(call_option,100.0_dp,105.0_dp,0.03_dp, &
                         1.5_dp,0.30_dp,0.02_dp)
   implied = implied_volatility(value%price,call_option,100.0_dp,105.0_dp, &
                                0.03_dp,1.5_dp,0.02_dp)
   write(*,'(a,f12.6)') 'Price:       ',value%price
   write(*,'(a,f12.6)') 'Delta:       ',value%delta
   write(*,'(a,f12.6)') 'Vega:        ',value%vega
   write(*,'(a,f12.6)') 'Implied vol: ',implied
end program example_black_scholes
