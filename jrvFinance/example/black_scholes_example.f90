! SPDX-License-Identifier: GPL-2.0-or-later
program black_scholes_example
  use jrvfinance, only: dp, black_scholes_result, gen_bs, gen_bs_implied
  implicit none
  type(black_scholes_result) :: result
  real(dp) :: implied

  result = gen_bs(100.0_dp, 100.0_dp, 0.10_dp, 0.20_dp, 1.0_dp)
  implied = gen_bs_implied(100.0_dp, 100.0_dp, 0.10_dp, result%call, 1.0_dp)
  print '(a,f10.4)', 'Call: ', result%call
  print '(a,f10.4)', 'Put: ', result%put
  print '(a,f10.6)', 'Delta: ', result%call_delta
  print '(a,f10.6)', 'Gamma: ', result%gamma
  print '(a,f10.6)', 'Implied volatility: ', implied
end program black_scholes_example
