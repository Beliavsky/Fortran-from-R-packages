program test_all_exports
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use bondanalyst
   implicit none

   integer :: failures, status
   real(dp) :: invalid_value
   real(dp), parameter :: times5(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   real(dp), parameter :: coupons5(5) = [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp]
   real(dp), parameter :: deficiency5(5) = [-2.0_dp, -2.0_dp, -2.0_dp, -2.0_dp, -2.0_dp]
   real(dp), parameter :: excess4(4) = [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
   real(dp), parameter :: times4(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   real(dp), parameter :: times6(6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
   real(dp), parameter :: times8(8) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
      5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
   real(dp), parameter :: spots3(3) = [0.02_dp, 0.03_dp, 0.04_dp]
   real(dp), parameter :: spotsz(3) = [0.0486_dp, 0.0495_dp, 0.0565_dp]
   real(dp), parameter :: curve8(8) = [0.10_dp, 0.12_dp, 0.15_dp, 0.18_dp, &
      0.20_dp, 0.22_dp, 0.24_dp, 0.30_dp]
   real(dp), parameter :: sacurve8(8) = [0.05_dp, 0.054_dp, 0.058_dp, 0.064_dp, &
      0.070_dp, 0.072_dp, 0.074_dp, 0.078_dp]
   real(dp), parameter :: fri3(3) = [1.0188_dp, 1.0188_dp*1.0277_dp, &
      1.0188_dp*1.0277_dp*1.0354_dp]

   failures = 0

   call check('discouponpmtsbond', discouponpmtsbond(coupons5, times5, 0.06_dp), &
      16.849_dp, 1.0e-12_dp)
   call check('dismaturityvalbond', dismaturityvalbond(100.0_dp, 5, 0.06_dp), &
      74.726_dp, 1.0e-12_dp)
   call check('bondpriceyearlycoupons', &
      bondpriceyearlycoupons(coupons5, times5, 100.0_dp, 5, 0.06_dp), &
      91.575_dp, 1.0e-12_dp)
   call check('pvcoupondeficiency', pvcoupondeficiency(deficiency5, times5, 0.06_dp), &
      -8.425_dp, 1.0e-12_dp)
   call check('bondpricedefcoupon', &
      bondpricedefcoupon(100.0_dp, deficiency5, times5, 0.06_dp), &
      91.575_dp, 1.0e-12_dp)
   call check('pvexcesscoupon', pvexcesscoupon(excess4, times4, 0.04_dp), &
      7.260_dp, 1.0e-12_dp)
   call check('bondpriceexcesscoupon', &
      bondpriceexcesscoupon(excess4, times4, 0.04_dp), 107.260_dp, 1.0e-12_dp)
   call check('pricingzerocouponbond', pricingzerocouponbond(100.0_dp, 10, 0.02_dp), &
      82.035_dp, 1.0e-12_dp)
   call check('ytmzerocouponbond', ytmzerocouponbond(100.0_dp, 60, 22.375_dp), &
      0.02527_dp, 1.0e-12_dp)
   call check('computingbondytmratefivedecimalplaces', &
      computingbondytmratefivedecimalplaces(5.0_dp, 100.0_dp, 105.0_dp, 4), &
      0.03634_dp, 1.0e-12_dp)
   call check('computingbondytmratesixdecimalplaces', &
      computingbondytmratesixdecimalplaces(5.0_dp, 100.0_dp, 105.0_dp, 4), &
      0.036344_dp, 1.0e-12_dp)
   call check('pricingsacpnbond', &
      pricingsacpnbond([4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], &
      times6, 100.0_dp, 3, 0.07_dp), 102.66_dp, 1.0e-12_dp)
   call check('pricingqtrlycpnbond', &
      pricingqtrlycpnbond([2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
      2.0_dp, 2.0_dp], times8, 100.0_dp, 2, 0.06_dp), 103.74_dp, 1.0e-12_dp)
   call check('pricingwithspots', &
      pricingwithspots([5.0_dp, 5.0_dp, 5.0_dp], spots3, &
      [1.0_dp, 2.0_dp, 3.0_dp], 100.0_dp, 3), 102.96_dp, 1.0e-12_dp)
   call check('pricingwithsptseq', &
      pricingwithsptseq([3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp], &
      [0.0408_dp, 0.0401_dp, 0.0370_dp, 0.0350_dp], times4, 100.0_dp, 4), &
      98.104_dp, 1.0e-12_dp)
   call check('aiactdtcon', aiactdtcon(4.375_dp, 0, 184, 43), &
      0.511209_dp, 1.0e-12_dp)
   call check('airoundeddaysconv', airoundeddaysconv(4.375_dp, 0, 12, 1, 180), &
      0.510417_dp, 1.0e-12_dp)
   call check('matrixmethod', &
      matrixmethod([4.5_dp, 4.5_dp, 4.5_dp, 4.5_dp], times4, &
      100.0_dp, 4, 0.02856_dp, 0.03449_dp), 104.991_dp, 1.0e-12_dp)
   call check('convertaprtodifferentperiodcity', &
      convertaprtodifferentperiodcity(0.0496_dp, 2, 4), 0.0493_dp, 1.0e-12_dp)
   call check('extracompensationforhigherrisk', &
      extracompensationforhigherrisk(0.10839_dp, 0.10630_dp), &
      20.9_dp, 1.0e-10_dp)
   call check('annualytmzcbforperiodicity', &
      annualytmzcbforperiodicity(100.0_dp, 5.0_dp, 80.0_dp, 2), &
      0.045130_dp, 1.0e-12_dp)
   call check('earzcbvariousperiodicity', &
      earzcbvariousperiodicity(100.0_dp, 5.0_dp, 80.0_dp, 2), &
      0.022565_dp, 1.0e-12_dp)
   call check('returnincomefrn', returnincomefrn(0.0125_dp, 0.0050_dp, 100.0_dp, 2), &
      0.875_dp, 1.0e-12_dp)
   call check('pricingfrn', &
      pricingfrn([0.875_dp, 0.875_dp, 0.875_dp, 0.875_dp], times4, &
      100.0_dp, 4, 0.00825_dp), 100.196_dp, 1.0e-12_dp)
   call check('periodicdiscratefrn', &
      periodicdiscratefrn(0.8125_dp, 100.0_dp, 98.0_dp, 4, 4), &
      0.009478_dp, 1.0e-12_dp)
   call check('discmarginfrn', discmarginfrn(0.0200_dp, 0.8125_dp, &
      100.0_dp, 98.0_dp, 4, 4), 0.017912_dp, 1.0e-12_dp)
   call check('computingytc', &
      computingytc(8.0_dp, 102.0_dp, 105.0_dp, 7, 4), 0.06975_dp, 1.0e-12_dp)
   call check('pricingtbill', pricingtbill(10000000.0_dp, 91, 360, 0.0225_dp), &
      9943125.0_dp, 1.0e-9_dp)
   call check('pricingmoneymarketinstrusingaor', &
      pricingmoneymarketinstrusingaor(10216000.0_dp, 180, 365, 0.0438_dp), &
      10000000.0_dp, 1.0e-9_dp)
   call check('fvmoneymarketinstrusingaor', &
      fvmoneymarketinstrusingaor(10000000.0_dp, 180, 365, 0.0438_dp), &
      10216000.0_dp, 1.0e-9_dp)
   call check('computingaormoneymarketinstr', &
      computingaormoneymarketinstr(10000000.0_dp, 10060829.0_dp, 45, 365), &
      0.04934_dp, 1.0e-12_dp)
   call check('pricingcommercialpaper', &
      pricingcommercialpaper(100.0_dp, 90, 360, 0.0576_dp), 98.56_dp, 1.0e-12_dp)
   call check('computingquoteddiscratemmi', &
      computingquoteddiscratemmi(9943125.0_dp, 10000000.0_dp, 91, 360), &
      0.02250_dp, 1.0e-12_dp)
   call check('fvmmiusingquoteddiscrate', &
      fvmmiusingquoteddiscrate(9943125.0_dp, 91, 360, 0.0225_dp), &
      10000000.0_dp, 1.0e-9_dp)
   call check('computingparrate', &
      computingparrate([0.05263_dp, 0.05616_dp, 0.06359_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp], 100.0_dp, 100.0_dp, 3), 6.306_dp, 1.0e-12_dp)
   call check('forwards', forwards(curve8, 1, 1, times8, 8), &
      0.1404_dp, 1.0e-12_dp)
   call check('saforwards', saforwards(sacurve8, 1, 1, times8, 8), &
      0.058008_dp, 1.0e-12_dp)
   call check('frpricing', frpricing([3.75_dp, 3.75_dp, 3.75_dp], fri3, &
      100.0_dp, 3), 102.965_dp, 1.0e-12_dp)
   call check('computinggspread', computinggspread(0.05932_dp, 0.03605_dp), &
      0.02327_dp, 1.0e-12_dp)
   call check('computingzspread', computingzspread(5.0_dp, 100.0_dp, &
      92.38_dp, 3, spotsz), 0.0234_dp, 1.0e-12_dp)
   call check('pricingwithzspread', &
      pricingwithzspread([5.0_dp, 5.0_dp, 5.0_dp], spotsz, &
      [1.0_dp, 2.0_dp, 3.0_dp], 100.0_dp, 3, 0.0234_dp), &
      92.383_dp, 1.0e-12_dp)
   call check('pricingwithgspread', &
      pricingwithgspread([6.0_dp, 6.0_dp], [1.0_dp, 2.0_dp], &
      100.0_dp, 2, 0.03605_dp, 0.02327_dp), 100.12_dp, 1.0e-12_dp)
   call check('macduration', macduration(10, 0.104_dp, 8.0_dp, 100.0_dp, 0, 0), &
      7.0029_dp, 1.0e-12_dp)
   call check('pvfullprice', pvfullprice(10, 0.104_dp, 8.0_dp, 100.0_dp, 0, 0), &
      85.5031_dp, 1.0e-12_dp)
   call check('macdurationonfp', &
      macdurationonfp(85.5031_dp, 10, 0.104_dp, 8.0_dp, 100.0_dp, 0, 0), &
      7.0029_dp, 1.0e-12_dp)
   call check('macdurationoncouponrate', &
      macdurationoncouponrate(0.08_dp, 10, 0.104_dp, 0.0_dp), &
      7.0029_dp, 1.0e-12_dp)
   call check('modifduration', &
      modifduration(10, 0.104_dp, 8.0_dp, 100.0_dp, 0, 0), &
      6.3432_dp, 1.0e-12_dp)
   call check('modifdurationusingmacduration', &
      modifdurationusingmacduration(7.0029_dp, 0.104_dp), 6.3432_dp, 1.0e-12_dp)
   call check('estimatedpercentchangepvfullprice', &
      estimatedpercentchangepvfullprice(6.126829_dp, 0.01_dp), &
      0.061268_dp, 1.0e-12_dp)
   call check('approxmodifduration', &
      approxmodifduration(99.956780_dp, 100.631781_dp, 101.250227_dp, 0.0005_dp), &
      6.187134079349011_dp, 1.0e-12_dp)
   call check('approxmacdurationusingapprmodifduration', &
      approxmacdurationusingapprmodifduration(13.466_dp, 0.0257_dp), &
      13.812_dp, 1.0e-12_dp)
   call check('effdurtncallablebond', &
      effdurtncallablebond(101.060489_dp, 99.050120_dp, 102.890738_dp, 0.0025_dp), &
      7.6006_dp, 1.0e-12_dp)
   call check('moneyduration', moneyduration(7.0029_dp, 0.104_dp, 85.5031_dp), &
      542.3638_dp, 1.0e-10_dp)
   call check('changepvfullbondprice', changepvfullbondprice(542.3638_dp, 0.01_dp), &
      -5.4236_dp, 1.0e-12_dp)
   call check('computingbondpvbp', computingbondpvbp(100.594327_dp, 100.765123_dp), &
      0.0854_dp, 1.0e-12_dp)

   call check('effectiveannualratezcb', &
      effectiveannualratezcb(100.0_dp, 5.0_dp, 80.0_dp), &
      (100.0_dp/80.0_dp)**0.2_dp-1.0_dp, 1.0e-13_dp)
   call check('conventionalpercentpricechange', &
      conventionalpercentpricechange(6.126829_dp, 0.01_dp), -0.06126829_dp, 1.0e-13_dp)

   invalid_value = discouponpmtsbond([1.0_dp], [1.0_dp, 2.0_dp], 0.05_dp, status)
   call check_nan('size mismatch value', invalid_value)
   call check_integer('size mismatch status', status, ba_size_mismatch)
   invalid_value = computingytc(8.0_dp, 100.0_dp, 105.0_dp, 5, 6, status)
   call check_nan('out of range value', invalid_value)
   call check_integer('out of range status', status, ba_out_of_range)

   if (failures /= 0) then
      error stop 'test_all_exports failed'
   end if
   print '(a)', 'test_all_exports: all checks passed'

contains

   subroutine check(name, actual, expected, tolerance)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance

      if (abs(actual-expected) > tolerance) then
         failures = failures+1
         print '(a,2(1x,es24.16))', trim(name)//' failed:', actual, expected
      end if
   end subroutine check


   subroutine check_nan(name, actual)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual

      if (.not. ieee_is_nan(actual)) then
         failures = failures+1
         print '(a,1x,es24.16)', trim(name)//' failed:', actual
      end if
   end subroutine check_nan


   subroutine check_integer(name, actual, expected)
      character(len=*), intent(in) :: name
      integer, intent(in) :: actual, expected

      if (actual /= expected) then
         failures = failures+1
         print '(a,2(1x,i0))', trim(name)//' failed:', actual, expected
      end if
   end subroutine check_integer

end program test_all_exports
