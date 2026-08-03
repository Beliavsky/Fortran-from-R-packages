program test_calibration_greeks
   use ragtop
   implicit none
   type(market_spec) :: market
   type(option_value) :: value
   type(instrument_spec) :: option
   type(greek_result) :: greeks
   real(dp) :: iv, jump_vol, bs_vol, american_price
   integer :: status

   value = black_scholes(put_option,100.0_dp,90.0_dp,0.03_dp,1.0_dp, &
                         0.45_dp,0.07_dp)
   iv = implied_volatility(value%price,put_option,100.0_dp,90.0_dp, &
                           0.03_dp,1.0_dp,0.07_dp,status=status)
   call assert_int(status,ragtop_ok,'Black-Scholes implied status')
   call assert_close(iv,0.45_dp,3.0e-8_dp,'Black-Scholes implied vol')

   market%short_rate = 0.0_dp
   market%default_intensity = 0.07_dp
   jump_vol = equivalent_jump_vola_to_bs(0.75_dp,20.0_dp,market,status)
   call assert_close(jump_vol,0.5649847_dp,5.0e-5_dp, &
                     'equivalent jump volatility')
   bs_vol = equivalent_bs_vola_to_jump(jump_vol,20.0_dp,market,status)
   call assert_close(bs_vol,0.75_dp,5.0e-7_dp,'equivalent BS volatility')

   market%short_rate = 0.02_dp
   market%default_intensity = 0.05_dp
   market%volatility = 0.55_dp
   american_price = american(put_option,33.0_dp,30.0_dp,0.77_dp,market,100)
   iv = american_implied_volatility(american_price,put_option,33.0_dp, &
                                    30.0_dp,0.77_dp,market,60,status=status)
   call assert_close(iv,0.55_dp,1.0e-3_dp,'American implied volatility')

   option = EuropeanOption(1.0_dp,100.0_dp,call_option)
   market%short_rate = 0.03_dp
   market%default_intensity = 0.0_dp
   market%volatility = 0.25_dp
   greeks = find_greeks(100.0_dp,120,option,market)
   call assert_true(greeks%delta > 0.0_dp .and. greeks%delta < 1.0_dp, &
                    'delta range')
   call assert_true(greeks%gamma > 0.0_dp,'gamma sign')
   call assert_true(greeks%vega > 0.0_dp,'vega sign')

   print '(a)', 'test_calibration_greeks: PASS'
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
   subroutine assert_int(actual, expected, label)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: label
      if (actual /= expected) error stop trim(label)
   end subroutine assert_int
end program test_calibration_greeks
