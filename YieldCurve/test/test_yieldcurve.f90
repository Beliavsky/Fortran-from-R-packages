program test_yieldcurve
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use yieldcurve, only : dp, yc_success, yc_invalid_argument, beta1_spot, beta2_spot, &
      beta1_forward, beta2_forward, factor_beta1, factor_beta2, ns_rates, &
      svensson_rates, nelson_siegel_fit, svensson_fit
   implicit none

   integer :: failures

   failures = 0
   call test_factors(failures)
   call test_ns_rates(failures)
   call test_svensson_rates(failures)
   call test_ns_fit(failures)
   call test_svensson_fit(failures)
   call test_matrix_api(failures)
   call test_missing_and_errors(failures)

   if (failures /= 0) then
      write(*, '(a, i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write(*, '(a)') 'All YieldCurve tests passed.'

contains

   subroutine test_factors(failures)
      integer, intent(inout) :: failures

      call assert_close('factor beta1', factor_beta1(0.5_dp, 2.0_dp), &
         0.6321205588285577_dp, 2.0e-14_dp, failures)
      call assert_close('factor beta2', factor_beta2(0.5_dp, 2.0_dp), &
         0.2642411176571153_dp, 2.0e-14_dp, failures)
      call assert_close('spot beta1', beta1_spot(2.0_dp, 4.0_dp), &
         0.7869386805747332_dp, 2.0e-14_dp, failures)
      call assert_close('spot beta2', beta2_spot(2.0_dp, 4.0_dp), &
         0.1804080208620997_dp, 2.0e-14_dp, failures)
      call assert_close('forward beta1', beta1_forward(2.0_dp, 4.0_dp), &
         0.6065306597126334_dp, 2.0e-14_dp, failures)
      call assert_close('forward beta2', beta2_forward(2.0_dp, 4.0_dp), &
         0.3032653298563167_dp, 2.0e-14_dp, failures)
      call assert_close('zero limit beta1', beta1_spot(0.0_dp, 2.0_dp), 1.0_dp, 1.0e-15_dp, failures)
      call assert_close('zero limit beta2', beta2_spot(0.0_dp, 2.0_dp), 0.0_dp, 1.0e-15_dp, failures)
   end subroutine test_factors

   subroutine test_ns_rates(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: maturity(8) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
         3.0_dp, 5.0_dp, 7.0_dp, 10.0_dp]
      real(dp), parameter :: coeff(4) = [4.0_dp, -1.5_dp, 2.0_dp, 0.55_dp]
      real(dp), parameter :: expected(8) = [ &
         2.7240791181978374_dp, 2.9179973476864838_dp, 3.2306914607113115_dp, &
         3.6374982491047136_dp, 3.8607335437825060_dp, 4.0423392119126380_dp, &
         4.0845470548585210_dp, 4.0823640233559390_dp]
      real(dp) :: curve(8)
      integer :: stat

      call ns_rates(coeff, maturity, curve, stat)
      call assert_true('NS rates status', stat == yc_success, failures)
      call assert_vector_close('NS rates values', curve, expected, 3.0e-13_dp, failures)
   end subroutine test_ns_rates

   subroutine test_svensson_rates(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: maturity(11) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
         3.0_dp, 5.0_dp, 7.0_dp, 10.0_dp, 15.0_dp, 20.0_dp, 30.0_dp]
      real(dp), parameter :: coeff(6) = [4.2_dp, -2.0_dp, 1.5_dp, 0.8_dp, 1.8_dp, 8.5_dp]
      real(dp), parameter :: expected_spot(11) = [ &
         2.4392184740618785_dp, 2.6498662670088910_dp, 2.9992737304964456_dp, &
         3.4849262382794066_dp, 3.7853191404284066_dp, 4.0984652226065940_dp, &
         4.2373421919609950_dp, 4.3281686041164030_dp, 4.3783699821298210_dp, &
         4.3865767070124955_dp, 4.3665624983665660_dp]
      real(dp), parameter :: expected_forward(11) = [ &
         2.6635156478259780_dp, 3.0450507142587580_dp, 3.6142925164716880_dp, &
         4.2390385782723140_dp, 4.4928242078177450_dp, 4.5960364481664840_dp, &
         4.5676053812382700_dp, 4.5147102012565800_dp, 4.4442704342479890_dp, &
         4.3792104016336800_dp, 4.2827933075849190_dp]
      real(dp) :: curve(11)
      integer :: stat

      call svensson_rates(coeff, maturity, curve, 'Spot', stat)
      call assert_true('Svensson spot status', stat == yc_success, failures)
      call assert_vector_close('Svensson spot values', curve, expected_spot, 5.0e-13_dp, failures)

      call svensson_rates(coeff, maturity, curve, stat=stat)
      call assert_true('Svensson forward status', stat == yc_success, failures)
      call assert_vector_close('Svensson forward values', curve, expected_forward, 5.0e-13_dp, failures)
   end subroutine test_svensson_rates

   subroutine test_ns_fit(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: maturity(8) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
         3.0_dp, 5.0_dp, 7.0_dp, 10.0_dp]
      real(dp), parameter :: source_coeff(4) = [4.0_dp, -1.5_dp, 2.0_dp, 0.55_dp]
      real(dp) :: rate(8), fitted(8), coeff(4)
      integer :: stat

      call ns_rates(source_coeff, maturity, rate, stat)
      call nelson_siegel_fit(rate, maturity, coeff, stat)
      call assert_true('NS fit status', stat == yc_success, failures)
      call ns_rates(coeff, maturity, fitted, stat)
      call assert_true('NS fitted rate status', stat == yc_success, failures)
      call assert_true('NS fit RMSE', sqrt(sum((fitted - rate)**2) / size(rate)) < 2.0e-3_dp, failures)
      call assert_true('NS beta0 constraint', coeff(1) > 0.0_dp .and. coeff(1) < 20.0_dp, failures)
   end subroutine test_ns_fit

   subroutine test_svensson_fit(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: maturity(11) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
         3.0_dp, 5.0_dp, 7.0_dp, 10.0_dp, 15.0_dp, 20.0_dp, 30.0_dp]
      real(dp), parameter :: source_coeff(6) = [4.2_dp, -2.0_dp, 1.5_dp, 0.8_dp, 1.8_dp, 8.5_dp]
      real(dp) :: rate(11), fitted(11), coeff(6)
      integer :: stat

      call svensson_rates(source_coeff, maturity, rate, 'spot', stat)
      call svensson_fit(rate, maturity, coeff, stat)
      call assert_true('Svensson fit status', stat == yc_success, failures)
      call svensson_rates(coeff, maturity, fitted, 'spot', stat)
      call assert_true('Svensson fitted rate status', stat == yc_success, failures)
      call assert_true('Svensson fit RMSE', sqrt(sum((fitted - rate)**2) / size(rate)) < 8.0e-3_dp, failures)
   end subroutine test_svensson_fit

   subroutine test_matrix_api(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: maturity(4) = [0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
      real(dp) :: coeff(2, 4), curves(2, 4), one(4)
      integer :: stat

      coeff(1, :) = [4.0_dp, -1.0_dp, 1.5_dp, 0.4_dp]
      coeff(2, :) = [3.0_dp, -0.5_dp, 0.8_dp, 0.7_dp]
      call ns_rates(coeff, maturity, curves, stat)
      call assert_true('matrix NS status', stat == yc_success, failures)
      call ns_rates(coeff(2, :), maturity, one, stat)
      call assert_vector_close('matrix NS row', curves(2, :), one, 1.0e-14_dp, failures)
   end subroutine test_matrix_api

   subroutine test_missing_and_errors(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: maturity(8) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
         3.0_dp, 5.0_dp, 7.0_dp, 10.0_dp]
      real(dp), parameter :: source_coeff(4) = [4.0_dp, -1.5_dp, 2.0_dp, 0.55_dp]
      real(dp) :: rate(8), coeff(4), curve(8), bad_maturity(8)
      integer :: stat

      call ns_rates(source_coeff, maturity, rate, stat)
      rate(3) = ieee_value(rate(3), ieee_quiet_nan)
      call nelson_siegel_fit(rate, maturity, coeff, stat)
      call assert_true('NS missing-value fit', stat == yc_success, failures)

      bad_maturity = maturity
      bad_maturity(4) = bad_maturity(3)
      call ns_rates(source_coeff, bad_maturity, curve, stat)
      call assert_true('invalid maturity status', stat == yc_invalid_argument, failures)
   end subroutine test_missing_and_errors

   subroutine assert_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tolerance
      integer, intent(inout) :: failures

      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         write(*, '(a, 2(1x, es24.16))') 'FAIL '//trim(name)//':', actual, expected
      end if
   end subroutine assert_close

   subroutine assert_vector_close(name, actual, expected, tolerance, failures)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      integer, intent(inout) :: failures
      real(dp) :: error

      error = maxval(abs(actual - expected))
      if (error > tolerance) then
         failures = failures + 1
         write(*, '(a, 1x, es24.16)') 'FAIL '//trim(name)//' max error:', error
      end if
   end subroutine assert_vector_close

   subroutine assert_true(name, condition, failures)
      character(len=*), intent(in) :: name
      logical, intent(in) :: condition
      integer, intent(inout) :: failures

      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL '//trim(name)
      end if
   end subroutine assert_true

end program test_yieldcurve
