program test_identities
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use bondanalyst
   implicit none

   integer :: failures, status
   real(dp) :: y, price, z, duration_value, modified_value, invalid_value
   real(dp), parameter :: spots(3) = [0.0486_dp, 0.0495_dp, 0.0565_dp]
   real(dp), parameter :: settlement_price = 100.9404_dp

   failures = 0

   y = computingbondytmratesixdecimalplaces(1.25_dp, 100.0_dp, 98.175677_dp, 6)
   price = level_price(1.25_dp, 100.0_dp, 6, y)
   call check('ytm repricing', price, 98.175677_dp, 4.0e-4_dp)

   z = computingzspread(5.0_dp, 100.0_dp, 92.38_dp, 3, spots)
   price = zspread_price(5.0_dp, 100.0_dp, spots, z)
   call check('z-spread repricing', price, 92.38_dp, 4.0e-3_dp)

   duration_value = macduration(16, 0.03_dp, 3.0_dp, 100.0_dp, 57, 180)
   modified_value = modifduration(16, 0.03_dp, 3.0_dp, 100.0_dp, 57, 180)
   call check('modified-duration identity', modified_value, &
      round_local(duration_value/1.03_dp, 4), 1.0e-12_dp)
   call check('full-price settlement example', &
      pvfullprice(16, 0.03_dp, 3.0_dp, 100.0_dp, 57, 180), &
      settlement_price, 1.0e-4_dp)
   call check('duration from full price', &
      macdurationonfp(settlement_price, 16, 0.03_dp, 3.0_dp, &
      100.0_dp, 57, 180), duration_value, 1.0e-4_dp)

   call check('forward compounding identity', &
      (1.0_dp+0.10_dp)*(1.0_dp+forwards([0.10_dp, 0.12_dp], 1, 1, &
      [1.0_dp, 2.0_dp], 2)), (1.0_dp+0.12_dp)**2, 6.0e-5_dp)

   invalid_value = pricingtbill(100.0_dp, 0, 360, 0.05_dp, status)
   call check_true('invalid day count value', ieee_is_nan(invalid_value))
   call check_integer('invalid day count status', status, ba_invalid_argument)

   invalid_value = computingbondytmratefivedecimalplaces(1.0_dp, 100.0_dp, &
      200.0_dp, 2, status)
   call check_true('no positive root value', ieee_is_nan(invalid_value))
   call check_integer('no positive root status', status, ba_no_root)

   invalid_value = pricingwithspots([5.0_dp, 5.0_dp], [0.02_dp], &
      [1.0_dp, 2.0_dp], 100.0_dp, 2, status)
   call check_true('spot size mismatch value', ieee_is_nan(invalid_value))
   call check_integer('spot size mismatch status', status, ba_size_mismatch)

   if (failures /= 0) error stop 'test_identities failed'
   print '(a)', 'test_identities: all checks passed'

contains

   pure function level_price(coupon, principal, periods, rate) result(value)
      real(dp), intent(in) :: coupon, principal, rate
      integer, intent(in) :: periods
      real(dp) :: value
      integer :: i

      value = 0.0_dp
      do i = 1, periods-1
         value = value+coupon/(1.0_dp+rate)**i
      end do
      value = value+(coupon+principal)/(1.0_dp+rate)**periods
   end function level_price


   pure function zspread_price(coupon, principal, curve, spread) result(value)
      real(dp), intent(in) :: coupon, principal, curve(:), spread
      real(dp) :: value
      integer :: i, n

      n = size(curve)
      value = 0.0_dp
      do i = 1, n-1
         value = value+coupon/(1.0_dp+curve(i)+spread)**i
      end do
      value = value+(coupon+principal)/(1.0_dp+curve(n)+spread)**n
   end function zspread_price


   pure function round_local(x, digits) result(value)
      real(dp), intent(in) :: x
      integer, intent(in) :: digits
      real(dp) :: value, scale

      scale = 10.0_dp**digits
      value = anint(x*scale)/scale
   end function round_local


   subroutine check(name, actual, expected, tolerance)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance

      if (abs(actual-expected) > tolerance) then
         failures = failures+1
         print '(a,2(1x,es24.16))', trim(name)//' failed:', actual, expected
      end if
   end subroutine check


   subroutine check_true(name, condition)
      character(len=*), intent(in) :: name
      logical, intent(in) :: condition

      if (.not. condition) then
         failures = failures+1
         print '(a)', trim(name)//' failed'
      end if
   end subroutine check_true


   subroutine check_integer(name, actual, expected)
      character(len=*), intent(in) :: name
      integer, intent(in) :: actual, expected

      if (actual /= expected) then
         failures = failures+1
         print '(a,2(1x,i0))', trim(name)//' failed:', actual, expected
      end if
   end subroutine check_integer

end program test_identities
