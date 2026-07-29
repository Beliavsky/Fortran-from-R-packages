program curve_and_duration
   use bondanalyst, only : dp, pricingwithspots, forwards, &
      approxmodifduration, computingparrate
   implicit none

   real(dp), parameter :: coupons(3) = [5.0_dp, 5.0_dp, 5.0_dp]
   real(dp), parameter :: spots(3) = [0.02_dp, 0.03_dp, 0.04_dp]
   real(dp), parameter :: times3(3) = [1.0_dp, 2.0_dp, 3.0_dp]
   real(dp), parameter :: curve(8) = [0.10_dp, 0.12_dp, 0.15_dp, 0.18_dp, &
      0.20_dp, 0.22_dp, 0.24_dp, 0.30_dp]
   real(dp), parameter :: times8(8) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
      5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]

   print '(a,f10.2)', 'spot-priced bond:       ', &
      pricingwithspots(coupons, spots, times3, 100.0_dp, 3)
   print '(a,f10.4)', '1y forward in one year:', forwards(curve, 1, 1, times8, 8)
   print '(a,f10.3)', 'three-year par rate %: ', &
      computingparrate(spots, times3, 100.0_dp, 100.0_dp, 3)
   print '(a,f10.6)', 'approx. mod. duration: ', &
      approxmodifduration(99.956780_dp, 100.631781_dp, 101.250227_dp, 0.0005_dp)
end program curve_and_duration
