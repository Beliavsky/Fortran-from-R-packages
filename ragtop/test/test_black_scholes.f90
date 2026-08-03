program test_black_scholes
   use ragtop
   implicit none
   type(option_value) :: value, direct
   type(market_spec) :: market
   type(dividend_schedule) :: dividends

   value = black_scholes(put_option,100.0_dp,90.0_dp,0.0_dp,1.0_dp,0.5_dp)
   call assert_close(value%price,14.160002955196_dp,1.0e-10_dp,'put price')
   call assert_close(value%delta,-0.322499382043183_dp,1.0e-10_dp,'put delta')
   call assert_close(value%vega,35.8771182787839_dp,1.0e-10_dp,'put vega')

   value = black_scholes(put_option,100.0_dp,90.0_dp,0.0_dp,1.0_dp, &
                         0.5_dp,0.07_dp)
   call assert_close(value%price,17.2748000186925_dp,1.0e-10_dp, &
                     'defaultable put price')

   dividends = make_dividend_schedule([0.5_dp],[0.0_dp],[20.0_dp])
   value = black_scholes(put_option,100.0_dp,90.0_dp,0.0_dp,1.0_dp, &
                         0.5_dp,0.07_dp,dividends=dividends)
   call assert_close(value%price,24.2940532694562_dp,1.0e-9_dp, &
                     'dividend put price')

   market%short_rate = 0.03_dp
   market%default_intensity = 0.07_dp
   market%volatility = 0.45_dp
   value = black_scholes_on_term_structures(put_option,100.0_dp,90.0_dp, &
                                             1.0_dp,market)
   direct = black_scholes(put_option,100.0_dp,90.0_dp,0.03_dp,1.0_dp, &
                           0.45_dp,0.07_dp)
   call assert_close(value%price,direct%price,1.0e-12_dp,'term structure')

   print '(a)', 'test_black_scholes: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
end program test_black_scholes
