! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1988-1990 Blake LeBaron
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program nonlinear_tests_example
  use fnonlinear, only : dp, rng_state, rng_seed, fill_normal, bds_test_result, &
    neural_test_result, runs_test_result, bds_test, white_neural_test, &
    terasvirta_neural_test, runs_test
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  real(dp) :: x(400), sd
  type(rng_state) :: rng
  type(bds_test_result) :: bds
  type(neural_test_result) :: white, tera
  type(runs_test_result) :: runs
  integer :: status

  call rng_seed(rng, 712367_int64)
  call fill_normal(rng, x)
  sd = sample_sd(x)
  call bds_test(x, 4, [0.5_dp * sd, sd, 1.5_dp * sd, 2.0_dp * sd], bds, status)
  call require(status == 0, "BDS failed")
  call white_neural_test(x, 2, 2, 10, 4.0_dp, white, status, rng=rng)
  call require(status == 0, "White test failed")
  call terasvirta_neural_test(x, 2, tera, status)
  call require(status == 0, "Terasvirta test failed")
  call runs_test(x, runs, status)
  call require(status == 0, "runs test failed")

  write(*, '(a, *(f10.4, 1x))') "BDS m=2..4, epsilon=sd: ", bds%statistic(:, 2)
  write(*, '(a, es14.6)') "White chi-square p-value: ", white%chi_square_p
  write(*, '(a, es14.6)') "Terasvirta chi-square p-value: ", tera%chi_square_p
  write(*, '(a, es14.6)') "Runs p-value: ", runs%p_value
contains
  pure real(dp) function sample_sd(values) result(value)
    real(dp), intent(in) :: values(:)
    real(dp) :: mean_value
    mean_value = sum(values) / real(size(values), dp)
    value = sqrt(sum((values - mean_value)**2) / real(size(values) - 1, dp))
  end function sample_sd

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) error stop message
  end subroutine require
end program nonlinear_tests_example
