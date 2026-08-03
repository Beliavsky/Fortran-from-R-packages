! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_core_hc
   use sandwich, only : dp, ols_model, fit_ols, meat, meat_hc, vcov_hc, vcov_opg, &
      SANDWICH_SUCCESS
   implicit none

   real(dp) :: x(6, 2), y(6)
   real(dp), allocatable :: m(:, :), mhc(:, :), cov(:, :), opg(:, :)
   type(ols_model) :: model
   integer :: status

   x(:, 1) = 1.0_dp
   x(:, 2) = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   y = [1.2_dp, 1.9_dp, 3.1_dp, 3.8_dp, 5.2_dp, 5.9_dp]

   call fit_ols(x, y, model, status)
   call assert_true(status == SANDWICH_SUCCESS, 'fit_ols status')
   call assert_close(model%coefficients(1), 3.02952380952381_dp, 1.0e-12_dp, 'intercept')
   call assert_close(model%coefficients(2), 0.974285714285714_dp, 1.0e-12_dp, 'slope')
   call assert_close(model%bread(1, 1), 1.08571428571429_dp, 1.0e-12_dp, 'bread 11')

   call meat(model%scores, m, status)
   call assert_true(status == SANDWICH_SUCCESS, 'meat status')
   call assert_close(m(1, 1), 0.0227936507936508_dp, 1.0e-12_dp, 'meat 11')
   call assert_close(m(1, 2), 0.0159682539682540_dp, 1.0e-12_dp, 'meat 12')
   call assert_close(m(2, 2), 0.0573312169312169_dp, 1.0e-12_dp, 'meat 22')

   call meat_hc(x, model%residuals, 'HC3', mhc, status, model%hat)
   call assert_true(status == SANDWICH_SUCCESS, 'HC3 meat status')
   call assert_close(mhc(1, 1), 0.048597138998001443_dp, 1.0e-11_dp, 'HC3 meat 11')
   call assert_close(mhc(1, 2), 0.020496788503528419_dp, 1.0e-11_dp, 'HC3 meat 12')
   call assert_close(mhc(2, 2), 0.14431641353396418_dp, 1.0e-11_dp, 'HC3 meat 22')

   call vcov_hc(x, model%residuals, model%bread, 'HC0', cov, status, model%hat)
   call assert_true(status == SANDWICH_SUCCESS, 'HC0 covariance status')
   call assert_close(cov(1, 1), 0.003768222006262827_dp, 1.0e-12_dp, 'HC0 vcov 11')
   call assert_close(cov(1, 2), -0.00019978490443796649_dp, 1.0e-12_dp, 'HC0 vcov 12')
   call assert_close(cov(2, 2), 0.00092201878846776811_dp, 1.0e-12_dp, 'HC0 vcov 22')

   call vcov_opg(model%scores, opg, status)
   call assert_true(status == SANDWICH_SUCCESS, 'OPG status')
   call assert_close(opg(1, 1), 9.08459704318817_dp, 1.0e-10_dp, 'OPG 11')
   call assert_close(opg(1, 2), -2.530299556678163_dp, 1.0e-10_dp, 'OPG 12')
   call assert_close(opg(2, 2), 3.6118391286181382_dp, 1.0e-10_dp, 'OPG 22')

   print '(a)', 'test_core_hc: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         print '(a,2(1x,es24.16))', 'FAIL: ' // trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_core_hc
