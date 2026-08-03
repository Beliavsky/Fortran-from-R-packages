program demo_fer
   use fer, only : dp, black_scholes_price, black_scholes_impvol, sabr_hagan_2002
   implicit none
   real(dp) :: price, iv
   price = black_scholes_price(100.0_dp,105.0_dp,1.2_dp,0.25_dp,1,exp(-0.04_dp*1.2_dp))
   iv = black_scholes_impvol(price,100.0_dp,105.0_dp,1.2_dp,1,exp(-0.04_dp*1.2_dp))
   print '(a,f12.6)', 'Black-Scholes call price: ', price
   print '(a,f12.6)', 'Recovered implied vol:   ', iv
   print '(a,f12.6)', 'SABR Hagan volatility:   ', &
      sabr_hagan_2002(1.0_dp,1.0_dp,10.0_dp,0.25_dp,0.3_dp,-0.8_dp,0.3_dp)
end program demo_fer
