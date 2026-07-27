! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
program test_inference
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchx_kinds, only : dp
   use garchx_math, only : set_random_seed, rmnorm, mean_value
   use garchx_inference, only : confidence_intervals, boundary_t_tests, boundary_wald_test
   implicit none
   integer, parameter :: n_draw = 30000
   integer :: status
   real(dp) :: par(3), vcov(3, 3), rmat(2, 3), restrictions(2), levels(3)
   real(dp) :: statistic
   real(dp), allocatable :: draws(:, :), ci(:, :), table(:, :), critical(:)
   real(dp) :: covariance(2, 2), means(2), sample_cov

   call set_random_seed(9031)
   covariance = reshape([1.0_dp, 0.35_dp, 0.35_dp, 2.0_dp], [2, 2])
   means = [0.5_dp, -0.2_dp]
   call rmnorm(n_draw, means, covariance, draws, status)
   call assert_true(status == 0, 'rmnorm status')
   call assert_true(abs(mean_value(draws(:, 1))-means(1)) < 0.025_dp, 'rmnorm mean 1')
   call assert_true(abs(mean_value(draws(:, 2))-means(2)) < 0.035_dp, 'rmnorm mean 2')
   sample_cov = sum((draws(:, 1)-mean_value(draws(:, 1)))* &
                    (draws(:, 2)-mean_value(draws(:, 2))))/real(n_draw-1, dp)
   call assert_true(abs(sample_cov-0.35_dp) < 0.035_dp, 'rmnorm covariance')

   par = [0.2_dp, 0.08_dp, 0.75_dp]
   vcov = 0.0_dp
   vcov(1, 1) = 0.01_dp
   vcov(2, 2) = 0.0016_dp
   vcov(3, 3) = 0.0064_dp
   vcov(2, 3) = -0.001_dp
   vcov(3, 2) = vcov(2, 3)
   call confidence_intervals(par, vcov, 500, 0.95_dp, ci, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ci)), 'confidence intervals')
   call assert_true(all(ci(:, 1) < par) .and. all(ci(:, 2) > par), 'confidence interval bounds')

   call boundary_t_tests(par, vcov, [2, 3], table, status)
   call assert_true(status == 0, 'boundary t test status')
   call assert_true(abs(table(1, 3)-2.0_dp) < 1.0e-12_dp, 'boundary t statistic')
   call assert_true(table(1, 4) > 0.02_dp .and. table(1, 4) < 0.03_dp, 'boundary t p-value')

   rmat = 0.0_dp
   rmat(1, 2) = 1.0_dp
   rmat(2, 3) = 1.0_dp
   restrictions = 0.0_dp
   levels = [0.10_dp, 0.05_dp, 0.01_dp]
   call boundary_wald_test(par, vcov, 500, rmat, restrictions, levels, 12000, &
                           statistic, critical, status)
   call assert_true(status == 0, 'boundary Wald status')
   call assert_true(statistic > 0.0_dp .and. all(critical > 0.0_dp), 'boundary Wald values')
   call assert_true(critical(1) < critical(2) .and. critical(2) < critical(3), &
                    'boundary Wald critical order')

   print '(a)', 'Multivariate-normal and boundary-inference tests passed.'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_inference
