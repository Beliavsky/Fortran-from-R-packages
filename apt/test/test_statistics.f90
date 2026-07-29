! SPDX-License-Identifier: GPL-2.0-or-later
program test_statistics
  use apt, only : dp, normal_cdf, student_t_cdf, f_cdf, chi_square_cdf, &
    regression_result, fit_ols
  use test_support, only : assert_close, assert_true
  implicit none
  real(dp) :: y(8), x(8,2)
  type(regression_result) :: fit
  integer :: i
  call assert_close(normal_cdf(1.25_dp), 0.8943502263331446_dp, 2e-13_dp, 'normal CDF')
  call assert_close(student_t_cdf(2.0_dp, 7.0_dp), 0.9571903357185120_dp, 2e-12_dp, 't CDF')
  call assert_close(f_cdf(3.4_dp, 4.0_dp, 15.0_dp), 0.9638992405479060_dp, 2e-12_dp, 'F CDF')
  call assert_close(chi_square_cdf(8.2_dp, 5.0_dp), 0.8544477570526422_dp, 2e-12_dp, 'chi-square CDF')
  do i = 1, 8
    x(i,1) = 1.0_dp
    x(i,2) = real(i,dp)
    y(i) = 1.25_dp + 0.75_dp*real(i,dp) + 0.1_dp*sin(real(i,dp))
  end do
  call fit_ols(y, x, fit)
  call assert_true(fit%status == 0, 'OLS status')
  call assert_close(fit%coefficients(1), 1.278342592492737_dp, 5e-12_dp, 'OLS intercept')
  call assert_close(fit%coefficients(2), 0.7479880099929075_dp, 5e-12_dp, 'OLS slope')
  print '(a)', 'test_statistics: PASS'
end program test_statistics
