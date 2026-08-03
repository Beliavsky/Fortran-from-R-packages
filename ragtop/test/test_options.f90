program test_options
   use ragtop
   implicit none
   type(market_spec) :: market
   type(instrument_spec) :: option
   real(dp) :: value

   market%short_rate = 0.06_dp
   market%volatility = 0.20_dp
   value = american(put_option,100.0_dp,110.0_dp,1.0_dp,market,200, &
                    std_devs_width=5.0_dp)
   call assert_close(value,11.6570723_dp,0.02_dp,'American control variate')

   market%short_rate = 0.02_dp
   market%volatility = 0.5_dp
   market%use_hazard_link = .true.
   market%hazard%base_intensity = 0.05_dp
   market%hazard%constant_fraction = 0.95_dp
   market%hazard%power = 1.0_dp
   market%hazard%reference_spot = 100.0_dp
   option = EuropeanOption(3.53_dp,200.0_dp,put_option)
   value = find_present_value(100.0_dp,250,option,market, &
                              std_devs_width=3.0_dp)
   call assert_close(value,109.2_dp,0.1_dp,'power hazard put')

   print '(a)', 'test_options: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
end program test_options
