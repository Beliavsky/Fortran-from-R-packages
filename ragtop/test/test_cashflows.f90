program test_cashflows
   use ragtop
   implicit none
   type(dividend_schedule) :: dividends
   type(cashflow_schedule) :: coupons
   type(market_spec) :: market
   real(dp) :: stock(6), values(6), shifted(6), expected(6)
   real(dp) :: prior, accelerated

   dividends = make_dividend_schedule([1.01_dp],[0.003_dp],[0.5_dp])
   stock = [11.0_dp,12.0_dp,13.0_dp,14.0_dp,15.0_dp,16.0_dp]
   call time_adjusted_dividends(dividends,1.2_dp,0.0_dp,0.0_dp, &
                                stock,12.0_dp,values)
   expected = [0.461333333333333_dp,0.503_dp,0.544666666666667_dp, &
               0.586333333333333_dp,0.628_dp,0.669666666666667_dp]
   call assert_vector(values,expected,1.0e-12_dp,'dividend values')

   values = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   call shift_for_dividends(values,stock,0.0_dp*stock,shifted)
   call assert_vector(shifted,values,1.0e-12_dp,'zero dividend shift')

   coupons = make_cashflow_schedule([0.5_dp,1.0_dp,1.5_dp], &
                                     [2.0_dp,2.0_dp,2.0_dp])
   market%short_rate = 0.04_dp
   prior = value_from_prior_coupons(1.0_dp,coupons,market)
   call assert_close(prior,2.0_dp*exp(0.04_dp*0.5_dp)+2.0_dp, &
                     1.0e-12_dp,'prior coupons')
   accelerated = accelerated_coupon_value(1.0_dp,coupons,market)
   call assert_close(accelerated,2.0_dp*exp(-0.04_dp*0.5_dp), &
                     1.0e-12_dp,'accelerated coupon')

   print '(a)', 'test_cashflows: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label
      if (maxval(abs(actual-expected)) > tolerance) then
         write(*,'(a,es24.14)') trim(label)//' max error: ', &
                                maxval(abs(actual-expected))
         error stop 1
      end if
   end subroutine assert_vector
end program test_cashflows
