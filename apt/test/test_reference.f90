! SPDX-License-Identifier: GPL-2.0-or-later
program test_reference
  use apt, only : dp, apt_tar, apt_mtar, ci_tar_fit_result, ci_tar_fit, &
    ecm_fit_result, ecm_symmetric_fit, ecm_asymmetric_fit
  use test_support, only : assert_close, assert_true, generate_prices
  implicit none
  real(dp) :: x(80), y(80)
  type(ci_tar_fit_result) :: fit
  type(ecm_fit_result) :: ecm
  integer :: i
  real(dp), parameter :: tar_coef(4) = [ &
    -0.110718939551750_dp, -0.056099890504295_dp, &
     1.417001156525542_dp, -0.814978255821078_dp]
  real(dp), parameter :: sym_x(6) = [ &
    0.015732170679111_dp, 1.846704543323182_dp, -0.944107052874196_dp, &
    0.021738237990791_dp, -0.043996157186036_dp, 0.022165445850104_dp]
  real(dp), parameter :: asy_y(11) = [ &
    0.022050756091324_dp, 0.534892818727375_dp, -0.167532375908309_dp, &
    0.649058077441346_dp, -0.007177283492377_dp, 1.462875975594321_dp, &
   -0.878471633923242_dp, 1.373491526104006_dp, -0.882238767199678_dp, &
    0.034559892763471_dp, -0.078938409361663_dp]

  call generate_prices(x, y)
  call ci_tar_fit(y, x, fit, apt_tar, 2, 0.0_dp)
  call assert_true(fit%status == 0, 'TAR fit status')
  call assert_close(fit%long_run%coefficients(1), 2.932096876699577_dp, 2e-11_dp, 'LR intercept')
  call assert_close(fit%long_run%coefficients(2), 1.326510319022199_dp, 2e-11_dp, 'LR slope')
  do i = 1, 4
    call assert_close(fit%threshold_regression%coefficients(i), tar_coef(i), 3e-10_dp, 'TAR coefficient')
  end do
  call assert_close(fit%sse, 0.009045626326643741_dp, 3e-10_dp, 'TAR SSE')
  call assert_close(fit%aic, -468.27797497136476_dp, 3e-10_dp, 'TAR AIC')
  call assert_close(fit%no_cointegration_test%f_statistic, 30.825611524213496_dp, 3e-9_dp, 'TAR no-CI F')
  call assert_close(fit%symmetry_test%f_statistic, 7.1819301871155785_dp, 3e-9_dp, 'TAR symmetry F')

  call ci_tar_fit(y, x, fit, apt_mtar, 1, 0.0_dp)
  call assert_true(fit%status == 0, 'MTAR fit status')
  call assert_close(fit%threshold_regression%coefficients(1), -0.315554829848687_dp, 3e-10_dp, 'MTAR coefficient 1')
  call assert_close(fit%threshold_regression%coefficients(2), -0.184360385530914_dp, 3e-10_dp, 'MTAR coefficient 2')
  call assert_close(fit%threshold_regression%coefficients(3), 0.794903437196556_dp, 3e-10_dp, 'MTAR coefficient 3')

  call ecm_symmetric_fit(y, x, ecm, 2)
  call assert_true(ecm%status == 0, 'symmetric ECM status')
  do i = 1, 6
    call assert_close(ecm%equation_x%coefficients(i), sym_x(i), 4e-10_dp, 'symmetric ECM x coefficient')
  end do

  call ecm_asymmetric_fit(y, x, ecm, 2, .true., apt_mtar, 0.0_dp)
  call assert_true(ecm%status == 0, 'asymmetric ECM status')
  do i = 1, 11
    call assert_close(ecm%equation_y%coefficients(i), asy_y(i), 5e-9_dp, 'asymmetric ECM y coefficient')
  end do
  print '(a)', 'test_reference: PASS'
end program test_reference
