! SPDX-License-Identifier: GPL-3.0-or-later
program test_curve_calibration
   use quant_bond_curves
   implicit none
   type(qbc_bond) :: bonds(3)
   type(qbc_calibration_result) :: fit, interpolation
   real(dp) :: terms(3), true_rates(3), prices(3), initial(3)
   integer :: i, st

   terms = [1.0_dp, 2.0_dp, 3.0_dp]
   true_rates = [0.03_dp, 0.04_dp, 0.05_dp]
   initial = 0.02_dp
   do i = 1, 3
      bonds(i)%analysis_date = make_date(2020,1,1)
      bonds(i)%maturity = make_date(2020+i,1,1)
      bonds(i)%frequency = 0
      bonds(i)%asset_type = qbc_asset_fixed
      bonds(i)%rate_type = qbc_rate_continuous
      prices(i) = exp(-true_rates(i)*terms(i))
   end do
   call bootstrap_curve(bonds, prices, terms, initial, fit, status=st)
   call assert_true(st == qbc_success, 'Rsolnp curve fit status')
   call assert_true(fit%objective < 1.0e-8_dp, 'curve fit objective')
   call assert_close(maxval(abs(fit%curve%rates-true_rates)), 0.0_dp, 3.0e-4_dp, 'recovered zero rates')

   call curve_calibration(terms, true_rates, [0.5_dp,1.5_dp,4.0_dp], interpolation, &
                          approximation=2, rate_type=qbc_rate_continuous, status=st)
   call assert_close(maxval(abs(interpolation%curve%rates-[0.03_dp,0.035_dp,0.05_dp])), 0.0_dp, 1.0e-14_dp, 'market curve interpolation')

   print '(a)', 'test_curve_calibration: PASS'
contains
   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) error stop label
   end subroutine
   subroutine assert_close(x, y, tol, label)
      real(dp), intent(in) :: x, y, tol
      character(len=*), intent(in) :: label
      if (abs(x-y) > tol) then
         write(*,*) label, x, y
         error stop 'assert_close'
      end if
   end subroutine
end program test_curve_calibration
