! SPDX-License-Identifier: MIT
program test_parametric
   use fixedincome
   implicit none
   type(interpolation_t) :: model
   type(fit_result_t) :: fit
   type(spot_rate_curve_t) :: curve
   type(term_t) :: terms
   real(dp) :: years(12), rates(12), query(3)
   real(dp), allocatable :: fitted(:)
   integer :: i, status

   years = [(real(i,dp)/2.0_dp, i=1,12)]
   rates = nelson_siegel(years, 0.06_dp, -0.025_dp, 0.035_dp, 0.8_dp)
   terms = term(years, 'years')
   curve = spotratecurve(rates, terms, 'continuous', 'actual/365', 'actual', status=status)
   model = interp_nelsonsiegel(0.05_dp, -0.01_dp, 0.01_dp, 0.4_dp)
   fit = fit_interpolation(model, curve, max_iterations=3000, tolerance=1.0e-12_dp)
   call check(fit%status == FI_OK .or. fit%status == FI_NO_CONVERGENCE, 'fit status')
   call check(fit%objective < 1.0e-8_dp, 'nelson siegel fit')
   call set_interpolation(curve, fit%model, status)
   fitted = interpolate(curve, curve%term_days, status)
   call check(sqrt(sum((fitted-rates)**2)/real(size(rates),dp)) < 3.0e-5_dp, 'fit interpolation')

   query = [0.0_dp, 1.0_dp, 10.0_dp]
   rates(1:3) = nelson_siegel_svensson(query, 0.05_dp, -0.02_dp, 0.03_dp, -0.01_dp, 0.7_dp, 2.0_dp)
   call check(all(abs(rates(1:3)) < 1.0_dp), 'nss evaluation')

   rates = nelson_siegel_svensson(years, 0.055_dp, -0.02_dp, 0.025_dp, -0.012_dp, 0.8_dp, 2.2_dp)
   curve = spotratecurve(rates, terms, 'continuous', 'actual/365', 'actual', status=status)
   model = interp_nelsonsiegelsvensson(0.05_dp, -0.01_dp, 0.01_dp, -0.005_dp, 0.5_dp, 1.5_dp)
   fit = fit_interpolation(model, curve, max_iterations=5000, tolerance=1.0e-11_dp)
   call check(fit%objective < 5.0e-8_dp, 'svensson fit')
   call check_close(nelson_siegel(0.0_dp,0.06_dp,-0.025_dp,0.035_dp,0.8_dp), &
      0.035_dp, 1.0e-12_dp, 'ns zero limit')

   print '(a)', 'test_parametric: PASS'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine check
   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual-expected) <= tolerance, label)
   end subroutine check_close
end program test_parametric
