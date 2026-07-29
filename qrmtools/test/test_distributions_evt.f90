! SPDX-License-Identifier: GPL-3.0-or-later
program test_distributions_evt
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use qrmtools, only : dp, dgev, pgev, qgev, dgpd, pgpd, qgpd, &
    dpar, ppar, qpar, dgpdtail, pgpdtail, qgpdtail, &
    fit_gpd_mom, fit_gpd_mle, fit_gev_mle, hill_estimator, &
    mean_excess_np, mean_excess_gpd, tail_estimator_gpd, &
    var_gpd, es_gpd, var_t01, es_t01, hill_result, fit_result
  implicit none

  real(dp), parameter :: tol = 2.0e-10_dp
  real(dp) :: parameters(2)
  real(dp), allocatable :: mean_excess(:,:)
  type(hill_result) :: hill
  type(fit_result) :: fit
  real(dp) :: xgpd(5)
  real(dp) :: xgev(10)

  call assert_close(dgev(3.0_dp, 0.2_dp, 1.0_dp, 2.0_dp), &
    0.11203386432543157_dp, tol)
  call assert_close(pgev(3.0_dp, 0.2_dp, 1.0_dp, 2.0_dp), &
    0.6690626526678188_dp, tol)
  call assert_close(qgev(0.9_dp, 0.2_dp, 1.0_dp, 2.0_dp), &
    6.684274065025338_dp, 5.0e-10_dp)
  call assert_close(pgev(qgev(0.73_dp, -0.15_dp, 0.3_dp, 1.2_dp), &
    -0.15_dp, 0.3_dp, 1.2_dp), 0.73_dp, tol)

  call assert_close(dgpd(2.0_dp, 0.2_dp, 1.5_dp), &
    0.1614115803251724_dp, tol)
  call assert_close(pgpd(2.0_dp, 0.2_dp, 1.5_dp), &
    0.6933179973821724_dp, tol)
  call assert_close(qgpd(0.95_dp, 0.2_dp, 1.5_dp), &
    6.1542315226956_dp, 5.0e-10_dp)
  call assert_close(var_gpd(0.95_dp, 0.2_dp, 1.5_dp), &
    6.1542315226956_dp, 5.0e-10_dp)
  call assert_close(es_gpd(0.95_dp, 0.2_dp, 1.5_dp), &
    9.567789403369499_dp, 8.0e-10_dp)

  call assert_close(dpar(2.0_dp, 3.0_dp, 2.0_dp), 0.09375_dp, tol)
  call assert_close(ppar(qpar(0.8_dp, 3.0_dp, 2.0_dp), 3.0_dp, 2.0_dp), &
    0.8_dp, tol)
  call assert_close(pgpdtail(qgpdtail(0.99_dp, 2.0_dp, 0.1_dp, &
    0.2_dp, 1.5_dp), 2.0_dp, 0.1_dp, 0.2_dp, 1.5_dp), 0.99_dp, tol)
  call assert_true(dgpdtail(3.0_dp, 2.0_dp, 0.1_dp, 0.2_dp, 1.5_dp) > 0.0_dp)

  xgpd = [0.2_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.5_dp]
  parameters = fit_gpd_mom(xgpd)
  call assert_close(parameters(1), -0.28955042527339014_dp, tol)
  call assert_close(parameters(2), 1.4700874848116647_dp, tol)

  fit = fit_gpd_mle(xgpd, estimate_covariance=.false., max_iterations=5000)
  call assert_true(fit%ok)
  call assert_true(fit%parameters(2) > 0.0_dp)
  call assert_true(ieee_is_finite(fit%log_likelihood))

  xgev = [-1.3_dp, -0.7_dp, -0.2_dp, 0.1_dp, 0.4_dp, &
    0.8_dp, 1.2_dp, 1.8_dp, 2.5_dp, 3.4_dp]
  fit = fit_gev_mle(xgev, estimate_covariance=.false., max_iterations=5000)
  call assert_true(fit%ok)
  call assert_true(fit%parameters(3) > 0.0_dp)

  hill = hill_estimator([10.0_dp, 8.0_dp, 7.0_dp, 5.0_dp, 4.0_dp, &
    3.0_dp, 2.5_dp, 2.0_dp, 1.5_dp, 1.0_dp], 2, 5)
  call assert_true(hill%ok)
  call assert_close(hill%tail_index(1), &
    1.0_dp / ((log(10.0_dp)+log(8.0_dp))/2.0_dp-log(8.0_dp)), tol)

  mean_excess = mean_excess_np([1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, &
    4.0_dp, 5.0_dp], 2)
  call assert_true(size(mean_excess,1) == 4)
  call assert_close(mean_excess(1,2), 2.5_dp, tol)
  call assert_close(mean_excess_gpd(2.0_dp, 0.2_dp, 1.5_dp), &
    2.375_dp, tol)
  call assert_close(tail_estimator_gpd(4.0_dp, 2.0_dp, 0.1_dp, &
    0.2_dp, 1.5_dp), 0.1_dp*(1.0_dp-pgpd(2.0_dp,0.2_dp,1.5_dp)), tol)

  call assert_close(var_t01(0.975_dp, 5.0_dp), 1.9911641278965473_dp, 2.0e-9_dp)
  call assert_close(es_t01(0.975_dp, 5.0_dp), 2.727802071641671_dp, 3.0e-9_dp)

  print '(a)', 'test_distributions_evt: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual-expected) > tolerance*max(1.0_dp,abs(expected))) then
      print *, 'mismatch:', actual, expected, abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_distributions_evt
