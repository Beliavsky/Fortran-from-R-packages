! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_nonparametric
   use sde
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none

   integer, parameter :: n = 201
   real(dp), parameter :: dt = 0.01_dp
   real(dp) :: x(n), data(n, 3), grid(5)
   real(dp), allocatable :: distances(:, :), operators(:, :, :), basis(:, :)
   type(kernel_estimate) :: drift_estimate, diffusion_estimate, density_estimate
   type(change_point_result) :: known_change, estimated_change
   real(dp) :: density_integral
   integer :: i, failures

   failures = 0
   call seed_rng(123456789_i64)

   x(1) = 0.0_dp
   do i = 2, n
      if (i <= 101) then
         x(i) = x(i-1)+0.10_dp*sqrt(dt)*random_normal()
      else
         x(i) = x(i-1)+0.50_dp*sqrt(dt)*random_normal()
      end if
   end do

   call kernel_drift(x, dt, drift_estimate, n_grid=64)
   call kernel_diffusion(x, dt, diffusion_estimate, n_grid=64)
   call kernel_density(x, density_estimate, n_grid=128)
   call assert_true(size(drift_estimate%x) == 64 .and. size(drift_estimate%y) == 64, &
      "kernel drift dimensions", failures)
   call assert_true(size(diffusion_estimate%x) == 64 .and. all(diffusion_estimate%y >= 0.0_dp), &
      "kernel diffusion values", failures)
   call assert_true(size(density_estimate%x) == 128 .and. all(density_estimate%y >= 0.0_dp), &
      "kernel density values", failures)
   density_integral = sum(0.5_dp*(density_estimate%y(1:127)+density_estimate%y(2:128))* &
      (density_estimate%x(2:128)-density_estimate%x(1:127)))
   call assert_true(abs(density_integral-1.0_dp) < 0.04_dp, &
      "kernel density approximately integrates to one", failures)

   call detect_change_point(x, dt, known_change, drift=zero_drift, diffusion=unit_diffusion, &
      theta=[real(dp) ::])
   call assert_true(known_change%index >= 65 .and. known_change%index <= 145, &
      "known-model change point location", failures)
   call assert_true(known_change%scale_after > 2.0_dp*known_change%scale_before, &
      "known-model change point scales", failures)

   call detect_change_point(x, dt, estimated_change)
   call assert_true(estimated_change%index >= 55 .and. estimated_change%index <= 155, &
      "nonparametric change point location", failures)
   call assert_true(estimated_change%scale_after > 1.5_dp*estimated_change%scale_before, &
      "nonparametric change point scales", failures)

   data(:, 1) = x
   data(:, 2) = x
   data(:, 3) = -0.6_dp*x
   do i = 2, n
      data(i, 3) = 0.25_dp*data(i-1, 3)+0.2_dp*random_normal()
   end do
   data(75, 2) = nan_dp()
   call markov_operator_distance(data, 12, distances, operators, spline_order=4)
   call assert_true(all(ieee_is_finite(distances)), "Markov distances finite", failures)
   call assert_true(maxval(abs(distances-transpose(distances))) < 1.0e-13_dp, &
      "Markov distance symmetry", failures)
   call assert_true(maxval(abs([distances(1, 1), distances(2, 2), distances(3, 3)])) < 1.0e-13_dp, &
      "Markov distance diagonal", failures)
   call assert_true(distances(1, 2) < 1.0e-3_dp, &
      "Markov distance tolerates interpolated duplicate", failures)
   call assert_true(distances(1, 3) > distances(1, 2), &
      "Markov distance separates different dynamics", failures)
   call assert_true(size(operators, 1) == 12 .and. size(operators, 3) == 3, &
      "Markov operator dimensions", failures)

   grid = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
   call bspline_basis_matrix(grid, 0.0_dp, 1.0_dp, 8, 3, basis)
   do i = 1, size(grid)
      call assert_close(sum(basis(i, :)), 1.0_dp, 2.0e-13_dp, &
         "B-spline partition of unity", failures)
      call assert_true(all(basis(i, :) >= -1.0e-14_dp), "B-spline nonnegative", failures)
   end do

   if (failures /= 0) then
      write(*, '(a, i0)') "test_nonparametric: failures = ", failures
      error stop 1
   end if
   write(*, '(a)') "test_nonparametric: PASS"

contains

   pure function zero_drift(t, state, theta) result(value)
      real(dp), intent(in) :: t, state, theta(:)
      real(dp) :: value
      value = 0.0_dp*(t+state+real(size(theta), dp))
   end function zero_drift

   pure function unit_diffusion(t, state, theta) result(value)
      real(dp), intent(in) :: t, state, theta(:)
      real(dp) :: value
      value = 1.0_dp+0.0_dp*(t+state+real(size(theta), dp))
   end function unit_diffusion

   subroutine assert_close(actual, target, tolerance, label, count)
      real(dp), intent(in) :: actual, target, tolerance
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. ieee_is_finite(actual) .or. abs(actual-target) > tolerance) then
         write(*, '(a, 2(1x, es23.15), a, es12.4)') "FAIL "//trim(label)//":", actual, target, &
            " tolerance=", tolerance
         count = count+1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label, count)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: count
      if (.not. condition) then
         write(*, '(a)') "FAIL "//trim(label)
         count = count+1
      end if
   end subroutine assert_true

end program test_nonparametric
