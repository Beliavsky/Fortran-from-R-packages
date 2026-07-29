! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program option_risk
  use dowd, only: dp, black_scholes_call_price, long_call_var, american_put_price_binomial, &
       american_put_var_binomial, american_put_es_binomial
  implicit none

  write(*,'(a,f12.6)') "European call: ", &
    black_scholes_call_price(100.0_dp,105.0_dp,0.03_dp,0.20_dp,0.5_dp)
  write(*,'(a,f12.6)') "Long call VaR: ", &
    long_call_var(100.0_dp,105.0_dp,0.03_dp,0.08_dp,0.20_dp,180.0_dp,0.99_dp,10.0_dp)
  write(*,'(a,f12.6)') "American put:  ", &
    american_put_price_binomial(100.0_dp,105.0_dp,0.03_dp,0.20_dp,180.0_dp,500)
  write(*,'(a,f12.6)') "American VaR:  ", &
    american_put_var_binomial(10000.0_dp,100.0_dp,105.0_dp,0.03_dp,0.20_dp,180.0_dp,100,0.99_dp,10.0_dp)
  write(*,'(a,f12.6)') "American ES:   ", &
    american_put_es_binomial(10000.0_dp,100.0_dp,105.0_dp,0.03_dp,0.20_dp,180.0_dp,100,0.99_dp,10.0_dp)
end program option_risk
