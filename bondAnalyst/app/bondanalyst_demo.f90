program bondanalyst_demo
   use bondanalyst, only : dp, bondpriceyearlycoupons, &
      computingbondytmratesixdecimalplaces, macduration, computingzspread
   implicit none

   real(dp), parameter :: coupons(5) = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
   real(dp), parameter :: times(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   real(dp), parameter :: spots(3) = [0.0486_dp, 0.0495_dp, 0.0565_dp]

   print '(a,f10.3)', 'annual-coupon bond price: ', &
      bondpriceyearlycoupons(coupons, times, 100.0_dp, 5, 0.06_dp)
   print '(a,f10.6)', 'yield per period:         ', &
      computingbondytmratesixdecimalplaces(5.0_dp, 100.0_dp, 105.0_dp, 4)
   print '(a,f10.4)', 'macaulay duration:        ', &
      macduration(10, 0.104_dp, 8.0_dp, 100.0_dp, 0, 0)
   print '(a,f10.4)', 'z-spread:                 ', &
      computingzspread(5.0_dp, 100.0_dp, 92.38_dp, 3, spots)
end program bondanalyst_demo
