! SPDX-License-Identifier: BSD-3-Clause
program test_shared_helpers
  use waveslim
  use waveslim_test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf
  use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
  implicit none

  real(dp) :: x(6)

  x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    ieee_value(0.0_dp,ieee_quiet_nan), ieee_value(0.0_dp,ieee_positive_inf)]
  call assert_close_scalar(mean_value(x),2.5_dp,1.0e-14_dp,'finite-only mean')
  call assert_close_scalar(variance_value(x),5.0_dp/3.0_dp,1.0e-14_dp,'finite-only sample variance')
  call assert_close_scalar(variance_value(x,.false.),1.25_dp,1.0e-14_dp,'finite-only population variance')
  call assert_close_scalar(median_value(x),2.5_dp,1.0e-14_dp,'finite-only median')
  call assert_close_scalar(mad_value(x),1.482602218505602_dp,1.0e-14_dp,'finite-only MAD')
  call assert_close_scalar(quantile_type7(x,0.25_dp),1.75_dp,1.0e-14_dp,'finite-only type-7 quantile')
  call assert_close_scalar(normal_cdf(0.0_dp),0.5_dp,1.0e-15_dp,'normal CDF')
  call assert_close_scalar(normal_quantile(0.975_dp),1.959963984540054_dp,1.0e-13_dp,'normal quantile')
  call assert_close_scalar(chi_square_cdf(2.0_dp,2.0_dp),1.0_dp-exp(-1.0_dp),1.0e-13_dp,'chi-square CDF')
  call assert_true(ieee_is_nan(median_value([ieee_value(0.0_dp,ieee_quiet_nan)])),'empty finite median')
  call assert_close_scalar(mean_value([ieee_value(0.0_dp,ieee_positive_inf)]),0.0_dp,0.0_dp,'empty finite mean')

  write(*,'(a)') 'test_shared_helpers: PASS'
end program test_shared_helpers
