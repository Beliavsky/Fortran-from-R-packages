program example_american_option
   use ragtop
   implicit none
   type(market_spec) :: market
   real(dp) :: price

   market%short_rate = 0.06_dp
   market%volatility = 0.20_dp
   price = american(put_option,100.0_dp,110.0_dp,1.0_dp,market,200, &
                    std_devs_width=5.0_dp)
   write(*,'(a,f12.6)') 'American put value: ',price
end program example_american_option
