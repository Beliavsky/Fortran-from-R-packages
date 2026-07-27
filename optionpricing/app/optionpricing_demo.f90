! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program optionpricing_demo
   use optionpricing, only : dp, european_result, greeks_result, bs_ec, &
      asian_call_app_lord, asian_call
   implicit none
   type(european_result) :: euro
   type(greeks_result) :: asian

   euro=bs_ec(0.25_dp,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   print '(a,f12.6)', 'European call price: ',euro%price
   print '(a,f12.6)', 'European call delta:',euro%delta
   print '(a,f12.6)', 'European call gamma:',euro%gamma
   print '(a,f12.6)', 'Lord Asian approximation: ', &
      asian_call_app_lord(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,.true.)

   asian=asian_call(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp, &
      'best','QMC',n=257,nout=12,a_gen=76,seed=4711)
   print '(a,3(f12.6,1x))', 'Asian QMC price/delta/gamma: ',asian%estimate
   print '(a,3(es12.4,1x))', '95% half-widths:             ',asian%error95
end program optionpricing_demo
