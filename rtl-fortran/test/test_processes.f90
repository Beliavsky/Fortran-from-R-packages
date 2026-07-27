! SPDX-License-Identifier: MIT
program test_processes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rtl, only: dp, path_result, ou_fit_result, multivariate_result
  use rtl, only: sim_gbm, sim_ou, sim_ou_time, sim_ou_jump, fit_ou, sim_multivariates
  implicit none

  type(path_result) :: path
  type(ou_fit_result) :: fit
  type(multivariate_result) :: multi
  real(dp) :: eps(4, 2), mu_time(2), mu_value(2), series(30)
  real(dp) :: prices(5, 2), a
  integer :: i

  eps = 0.0_dp
  path = sim_gbm(2, 10.0_dp, 0.05_dp, 0.2_dp, 1.0_dp, 0.25_dp, eps)
  call assert_true(path%status%ok)
  do i = 0, 4
    call assert_close(path%values(i, 1), 10.0_dp * exp(0.03_dp * 0.25_dp * real(i, dp)), 1.0e-12_dp)
  end do

  path = sim_ou(2, 5.0_dp, 7.0_dp, 0.5_dp, 0.2_dp, 1.0_dp, 0.25_dp, eps)
  call assert_close(path%values(1, 1), 5.25_dp, 1.0e-14_dp)
  call assert_close(path%values(2, 1), 5.46875_dp, 1.0e-14_dp)

  mu_time = [0.0_dp, 1.0_dp]
  mu_value = [2.0_dp, 4.0_dp]
  path = sim_ou_time(2, 0.0_dp, mu_time, mu_value, 1.0_dp, 0.0_dp, 1.0_dp, 0.25_dp, eps)
  call assert_close(path%values(1, 1), 0.625_dp, 1.0e-14_dp)

  a = exp(-0.5_dp / 12.0_dp)
  series(1) = 8.0_dp
  do i = 2, size(series)
    series(i) = 5.0_dp + a * (series(i - 1) - 5.0_dp)
  end do
  fit = fit_ou(series, 1.0_dp / 12.0_dp)
  call assert_true(fit%status%ok)
  call assert_close(fit%theta, 0.5_dp, 2.0e-10_dp)
  call assert_close(fit%mu, 5.0_dp, 2.0e-10_dp)
  call assert_close(fit%sigma, 0.0_dp, 2.0e-7_dp)

  prices(:, 1) = [10.0_dp, 11.0_dp, 13.0_dp, 14.0_dp, 16.0_dp]
  prices(:, 2) = [20.0_dp, 19.0_dp, 21.0_dp, 20.0_dp, 22.0_dp]
  multi = sim_multivariates(20, prices, s0=[1.0_dp, 2.0_dp], seed=123)
  call assert_true(multi%status%ok)
  call assert_close(multi%mean(1), 1.5_dp, 1.0e-14_dp)
  call assert_close(multi%mean(2), 0.5_dp, 1.0e-14_dp)
  call assert_close(multi%covariance(1, 2), multi%covariance(2, 1), 0.0_dp)
  call assert_true(all(ieee_is_finite(multi%simulations(:, 1))))

  path = sim_ou_jump(2, 5.0_dp, 5.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, &
    0.1_dp, 1.0_dp, 0.25_dp, seed=3)
  call assert_true(path%status%ok)
  call assert_close(path%values(4, 1), 5.0_dp, 1.0e-14_dp)

  print '(a)', 'test_processes: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print '(a,3es24.15)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_processes
