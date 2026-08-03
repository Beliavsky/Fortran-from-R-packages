! SPDX-License-Identifier: GPL-2.0-or-later
program demo_jrvfinance
  use jrvfinance, only: dp, date_t, parse_date, npv, irr, annuity_pv, &
    annuity_instalment, bond_price, bond_yield, black_scholes_result, &
    gen_bs, gen_bs_implied, equiv_rate, edate, date_string
  implicit none
  type(date_t) :: settle, mature
  type(black_scholes_result) :: bs
  real(dp) :: price

  print '(a,f12.4)', 'NPV: ', npv([100.0_dp,250.0_dp,300.0_dp],0.05_dp)
  print '(a,f12.6)', 'IRR: ', irr([-600.0_dp,300.0_dp,400.0_dp])
  print '(a,f12.4)', 'Annuity PV: ', annuity_pv(0.10_dp,15.0_dp)
  print '(a,f12.4)', 'Loan payment: ', annuity_instalment(0.09_dp,8.0_dp,10000.0_dp)

  settle=parse_date('2012-04-15');mature=parse_date('2022-01-01')
  price=bond_price(settle,mature,0.08_dp,0.088843_dp)
  print '(a,f12.4)', 'Bond price: ', price
  print '(a,f12.6)', 'Recovered yield: ', bond_yield(settle,mature,0.08_dp,price)

  bs=gen_bs(100.0_dp,100.0_dp,0.10_dp,0.20_dp,1.0_dp)
  print '(a,f12.4)', 'Black-Scholes call: ', bs%call
  print '(a,f12.6)', 'Implied volatility: ', &
    gen_bs_implied(100.0_dp,100.0_dp,0.10_dp,bs%call,1.0_dp)
  print '(a,f12.6)', '10% monthly -> semiannual: ', equiv_rate(0.10_dp,12.0_dp,2.0_dp)
  print '(a,a)', '2007-02-28 plus 4 months: ', &
    date_string(edate(parse_date('2007-02-28'),4))
end program demo_jrvfinance
