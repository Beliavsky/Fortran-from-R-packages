! SPDX-License-Identifier: GPL-2.0-or-later
program test_ecm_tests
  use apt, only : dp, apt_mtar, ecm_fit_result, ecm_asymmetry_test_result, &
    ecm_diagnostics_result, ecm_asymmetric_fit, ecm_asymmetry_tests, &
    ecm_diagnostics
  use test_support, only : assert_true, assert_close, generate_prices
  implicit none
  real(dp) :: x(80), y(80)
  type(ecm_fit_result) :: fit
  type(ecm_asymmetry_test_result) :: tests
  type(ecm_diagnostics_result) :: diag
  call generate_prices(x, y)
  call ecm_asymmetric_fit(y, x, fit, 2, .true., apt_mtar, 0.0_dp)
  call ecm_asymmetry_tests(fit, tests)
  call assert_true(tests%status == 0, 'asymmetry test status')
  call assert_true(tests%ntests == 9, 'asymmetry test count')
  call assert_close(tests%tests(1)%equation_x%f_statistic, &
    3.112284627685974_dp, 5e-8_dp, 'H1 x F')
  call assert_close(tests%tests(1)%equation_y%f_statistic, &
    18.997638530884746_dp, 5e-8_dp, 'H1 y F')
  call assert_close(tests%tests(2)%equation_x%f_statistic, &
    1605.4453642894446_dp, 5e-7_dp, 'H2 x-to-x F')
  call assert_close(tests%tests(3)%equation_y%f_statistic, &
    340.00042691845385_dp, 5e-7_dp, 'H2 y-to-y F')

  call ecm_diagnostics(fit, diag)
  call assert_true(diag%status == 0, 'diagnostic status')
  call assert_close(diag%equation_x%durbin_watson, 0.16613337236921194_dp, &
    5e-9_dp, 'DW x')
  call assert_close(diag%equation_y%durbin_watson, 0.8579032828778721_dp, &
    5e-9_dp, 'DW y')
  call assert_close(diag%equation_x%ljung_box_4_p, 1.585122081373518e-40_dp, &
    2e-12_dp, 'LB4 x p')
  print '(a)', 'test_ecm_tests: PASS'
end program test_ecm_tests
