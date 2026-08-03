program test_bonds
   use ragtop
   implicit none
   type(market_spec) :: market
   type(cashflow_schedule) :: coupons
   type(exercise_schedule) :: calls, puts
   type(instrument_spec) :: zero, coupon, callable, convertible
   real(dp) :: p_zero, p_coupon, p_callable, p_convertible

   market%short_rate = 0.06_dp
   market%volatility = 0.25_dp
   zero = ZeroCouponBond(2.0_dp,100.0_dp,0.4_dp)
   p_zero = find_present_value(100.0_dp,100,zero,market)
   call assert_close(p_zero,100.0_dp*exp(-0.12_dp),1.0e-8_dp,'zero bond')

   coupons = make_cashflow_schedule([0.5_dp,1.0_dp,1.5_dp,2.0_dp], &
                                     [2.0_dp,2.0_dp,2.0_dp,2.0_dp])
   coupon = CouponBond(2.0_dp,100.0_dp,coupons,0.4_dp)
   p_coupon = find_present_value(100.0_dp,100,coupon,market)
   call assert_true(p_coupon > p_zero,'coupon bond premium')

   calls = make_exercise_schedule([1.0_dp],[98.0_dp])
   puts = make_exercise_schedule([1.0_dp],[92.0_dp])
   callable = CallableBond(2.0_dp,100.0_dp,coupons,calls,puts,0.4_dp)
   p_callable = find_present_value(100.0_dp,100,callable,market)
   call assert_true(p_callable <= p_coupon+1.0e-8_dp,'callable cap')

   convertible = ConvertibleBond(2.0_dp,100.0_dp,1.0_dp,coupons,0.4_dp)
   p_convertible = find_present_value(100.0_dp,100,convertible,market)
   call assert_true(p_convertible >= p_coupon,'convertible over straight bond')
   call assert_true(p_convertible >= 100.0_dp,'conversion floor')
   call assert_close(p_convertible,117.88360868723908_dp,1.0e-8_dp, &
                     'convertible reference value')

   print '(a)', 'test_bonds: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) error stop trim(label)
   end subroutine assert_true
end program test_bonds
