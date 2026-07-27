! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1988-1990 Blake LeBaron
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program test_tests
  use fnonlinear, only : dp, rng_state, rng_seed, logistic_sim, bds_test_result, &
    neural_test_result, runs_test_result, generic_test_result, bds_test, &
    white_neural_test, terasvirta_neural_test, runs_test, ts_test
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  real(dp), parameter :: sample(12) = [0.1_dp, 0.4_dp, 0.2_dp, 0.8_dp, &
    0.3_dp, 0.7_dp, 0.5_dp, 0.9_dp, 0.6_dp, 0.05_dp, 0.45_dp, 0.25_dp]
  real(dp), allocatable :: x(:), gamma(:, :)
  type(bds_test_result) :: bds
  type(neural_test_result) :: white, tera
  type(runs_test_result) :: runs
  type(generic_test_result) :: generic
  type(rng_state) :: rng
  integer :: status

  call bds_test(sample, 3, [0.2_dp, 0.3_dp], bds, status)
  call check(status == 0, "BDS status")
  call check(maxval(abs(bds%correlation(:, 1) - [1.0_dp / 3.0_dp, &
    0.08888888888888889_dp, 0.04444444444444444_dp])) < 1.0e-14_dp, &
    "BDS correlations")
  call check(abs(bds%k_value(1) - 0.09722222222222222_dp) < 1.0e-14_dp, &
    "BDS k")
  call check(maxval(abs(bds%statistic(:, 1) - [-2.529822128134697_dp, &
    1.1457837984630885_dp])) < 1.0e-11_dp, "BDS statistics epsilon one")
  call check(maxval(abs(bds%statistic(:, 2) - [8.184718649847499_dp, &
    15.851039450217622_dp])) < 1.0e-10_dp, "BDS statistics epsilon two")
  call check(all(bds%p_value >= 0.0_dp .and. bds%p_value <= 1.0_dp), "BDS p-values")

  call logistic_sim(300, 100, 3.9_dp, x, status, start=0.123_dp)
  call check(status == 0, "logistic test series")
  call rng_seed(rng, 12345_int64)
  call white_neural_test(x, 2, 2, 10, 4.0_dp, white, status, rng=rng, gamma_out=gamma)
  call check(status == 0, "White test status")
  call check(abs(gamma(1, 1) + 1.61350258_dp) < 2.0e-8_dp, "White deterministic weights")
  call check(abs(white%ssr_null - 17.796392908951123_dp) < 2.0e-10_dp, "White null SSR")
  call check(abs(white%ssr_alternative - 1.3634229324198304_dp) < 2.0e-9_dp, &
    "White alternative SSR")
  call check(abs(white%chi_square - 770.699217410901_dp) < 2.0e-7_dp, &
    "White chi-square")
  call check(abs(white%f_statistic - 1783.8042024202482_dp) < 2.0e-6_dp, &
    "White F statistic")
  call check(white%chi_square_p < 1.0e-100_dp .and. white%f_p < 1.0e-100_dp, &
    "White p-values")

  call terasvirta_neural_test(x, 2, tera, status)
  call check(status == 0, "Terasvirta status")
  call check(tera%numerator_df == 7 .and. tera%denominator_df == 291, &
    "Terasvirta degrees of freedom")
  call check(tera%chi_square > 10000.0_dp .and. tera%f_statistic > 1.0e20_dp, &
    "Terasvirta nonlinear detection")
  call check(tera%chi_square_p < 1.0e-100_dp .and. tera%f_p < 1.0e-100_dp, &
    "Terasvirta p-values")

  call runs_test([-1.0_dp, -2.0_dp, 0.0_dp, 3.0_dp, 4.0_dp, -1.0_dp, 2.0_dp], &
    runs, status)
  call check(status == 0 .and. runs%runs == 4, "runs count")
  call check(runs%negative_count == 3 .and. runs%positive_count == 3, "runs counts")
  call check(abs(runs%statistic) < 1.0e-14_dp .and. abs(runs%p_value - 1.0_dp) < 1.0e-14_dp, &
    "runs statistic")

  call ts_test(sample, "bds", generic, status, max_dimension=3, epsilon=[0.2_dp])
  call check(status == 0 .and. size(generic%statistic) == 2, "generic BDS dispatcher")
  call ts_test(x, "tnn", generic, status, lag=2)
  call check(status == 0 .and. size(generic%statistic) == 2, "generic TNN dispatcher")

  print '(a)', "BDS, neural-network, dispatcher, and runs tests passed."
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // message
      error stop 1
    end if
  end subroutine check
end program test_tests
