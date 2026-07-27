! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program demo_fnonlinear
  use fnonlinear, only : dp, rng_state, rng_seed, logistic_sim, henon_sim, &
    mutual_information_curve, false_nearest_neighbors, bds_test_result, &
    neural_test_result, runs_test_result, bds_test, terasvirta_neural_test, runs_test
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  real(dp), allocatable :: x(:), xy(:, :), ami(:), fractions(:)
  integer, allocatable :: totals(:)
  type(bds_test_result) :: bds
  type(neural_test_result) :: neural
  type(runs_test_result) :: runs
  type(rng_state) :: rng
  real(dp) :: sd
  integer :: status

  call rng_seed(rng, 20260725_int64)
  call logistic_sim(500, 200, 3.9_dp, x, status, start=0.123_dp)
  call require(status == 0, "logistic simulation failed")
  call henon_sim(500, 200, 1.4_dp, 0.3_dp, xy, status, start=[0.1_dp, 0.2_dp])
  call require(status == 0, "Henon simulation failed")
  sd = sample_sd(x)
  call mutual_information_curve(x, 16, 10, ami, status)
  call require(status == 0, "mutual information failed")
  call false_nearest_neighbors(x, 5, 1, 5, 10.0_dp, 0.2_dp * sd, fractions, totals, &
    status, search_method="auto")
  call require(status == 0, "false-nearest calculation failed")
  call bds_test(x, 3, [0.5_dp * sd, sd], bds, status)
  call require(status == 0, "BDS test failed")
  call terasvirta_neural_test(x, 2, neural, status)
  call require(status == 0, "Terasvirta test failed")
  call runs_test(x - sum(x) / real(size(x), dp), runs, status)
  call require(status == 0, "runs test failed")

  write(*, '(a)') "fNonlinear modern Fortran demonstration"
  write(*, '(a, es14.6)') "logistic mean: ", sum(x) / real(size(x), dp)
  write(*, '(a, es14.6)') "logistic standard deviation: ", sd
  write(*, '(a, es14.6)') "Henon x mean: ", sum(xy(:, 1)) / real(size(xy, 1), dp)
  write(*, '(a, es14.6)') "mutual information lag 1: ", ami(1)
  write(*, '(a, *(f8.4, 1x))') "false-neighbor fractions: ", fractions
  write(*, '(a, *(f10.4, 1x))') "BDS statistics at epsilon=sd: ", bds%statistic(:, 2)
  write(*, '(a, es14.6)') "Terasvirta chi-square p-value: ", neural%chi_square_p
  write(*, '(a, es14.6)') "runs-test p-value: ", runs%p_value
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
end program demo_fnonlinear
