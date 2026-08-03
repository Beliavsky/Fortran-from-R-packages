! SPDX-License-Identifier: GPL-2.0-or-later
program cashflows_annuity
  use jrvfinance, only: dp, npv, irr, annuity_pv, annuity_instalment, annuity_rate
  implicit none
  real(dp), parameter :: cf(3) = [-600.0_dp, 300.0_dp, 400.0_dp]
  real(dp) :: rate

  rate = irr(cf)
  print '(a,f12.6)', 'IRR: ', rate
  print '(a,f12.4)', 'NPV at IRR: ', npv(cf, rate, cf_t=[0.0_dp,1.0_dp,2.0_dp])
  print '(a,f12.4)', 'PV of 15-year annuity: ', annuity_pv(0.10_dp, 15.0_dp)
  print '(a,f12.4)', 'Payment on 8-year loan: ', &
    annuity_instalment(0.09_dp, 8.0_dp, 10000.0_dp)
  print '(a,f12.6)', 'Mortgage rate: ', annuity_rate(360.0_dp, 450.0_dp, &
    50000.0_dp, cf_freq=12.0_dp, comp_freq=2.0_dp)
end program cashflows_annuity
