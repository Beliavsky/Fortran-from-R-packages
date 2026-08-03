program test_term_structures
   use ragtop
   implicit none
   type(market_spec) :: market
   integer :: status
   real(dp) :: value

   call initialize_volatility_curve(market%vols,[0.1_dp,2.0_dp,3.0_dp], &
                                    [0.2_dp,0.5_dp,1.2_dp],status)
   call assert_int(status,ragtop_ok,'vol curve status')
   market%use_vol_curve = .true.
   value = sqrt(cumulative_variance(market,0.1_dp,0.0_dp)/0.1_dp)
   call assert_close(value,0.2_dp,1.0e-12_dp,'short volatility')
   value = sqrt(cumulative_variance(market,100.0_dp,0.0_dp)/100.0_dp)
   call assert_close(value,1.93613_dp,1.0e-5_dp,'long volatility')

   call initialize_discount_curve(market%rates,[1.0_dp,5.0_dp,10.0_dp], &
                                  [0.01_dp,0.02_dp,0.03_dp],status)
   call assert_int(status,ragtop_ok,'rate curve status')
   market%use_rate_curve = .true.
   call assert_close(discount_factor(market,1.0_dp,0.0_dp),exp(-0.01_dp), &
                     1.0e-12_dp,'one-year discount')
   call assert_close(discount_factor(market,5.0_dp,1.0_dp), &
                     exp(-0.10_dp)/exp(-0.01_dp),1.0e-12_dp, &
                     'forward discount')

   call check_discount_factor(market,[0.1_dp,1.0_dp,5.0_dp,10.0_dp],status)
   call assert_int(status,ragtop_ok,'discount validation')
   call check_variance_cumulation(market,[0.1_dp,1.0_dp,5.0_dp],status)
   call assert_int(status,ragtop_ok,'variance validation')

   market%use_hazard_link = .true.
   market%hazard%base_intensity = 0.05_dp
   market%hazard%constant_fraction = 0.95_dp
   market%hazard%power = 1.0_dp
   market%hazard%reference_spot = 100.0_dp
   call assert_close(default_intensity(market,0.0_dp,50.0_dp),0.0525_dp, &
                     1.0e-12_dp,'linked hazard')

   print '(a)', 'test_term_structures: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_int(actual, expected, label)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: label
      if (actual /= expected) then
         write(*,'(a,2i8)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_int
end program test_term_structures
