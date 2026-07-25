! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
program test_core
   use garchx_kinds, only : dp
   use garchx_math, only : normal_cdf, normal_quantile, student_t_cdf, student_t_quantile
   use garchx_utilities, only : glag, gdiff
   use garchx_model, only : garchx_spec, make_garchx_spec, garchx_filter
   use garchx_model, only : garchx_filter_derivatives, garchx_objective_value
   use garchx_model, only : garchx_simulate
   implicit none
   type(garchx_spec) :: spec
   integer :: status, j
   real(dp) :: objective, objective_manual, h
   real(dp), allocatable :: lagged(:), differenced(:), sigma2(:), deriv(:, :)
   real(dp), allocatable :: matrix_lag(:, :), matrix_diff(:, :)
   real(dp), allocatable :: sigma_plus(:), sigma_minus(:), ysim(:), ssim(:), zout(:)
   real(dp) :: x(8), y(8), pars(7), expected(8), z(8)
   real(dp) :: pplus(7), pminus(7), fd(8)
   real(dp) :: xreg(8, 1), matrix_data(4, 2), yzero(8)

   y = [1.0_dp, -2.0_dp, 0.5_dp, -1.0_dp, 1.5_dp, -0.75_dp, 0.25_dp, -1.25_dp]
   x = [(0.1_dp*real(j, dp), j=1, 8)]
   xreg(:, 1) = x
   pars = [0.2_dp, 0.1_dp, 0.03_dp, 0.5_dp, 0.1_dp, 0.07_dp, 0.04_dp]
   expected = [1.5_dp, 1.5_dp, 1.832078125_dp, 1.3370390625_dp, &
               1.36172734375_dp, 1.271067578125_dp, 1.1253315234375_dp, &
               0.99552251953125_dp]
   call make_garchx_spec(spec, arch_lags=[1, 3], garch_lags=[1, 2], &
                         asym_lags=[1], xreg_count=1)
   call garchx_filter(y, spec, pars, sigma2, status, xreg, 1.5_dp)
   call assert_true(status == 0, 'filter status')
   call assert_close(maxval(abs(sigma2-expected)), 0.0_dp, 2.0e-12_dp, 'filter values')

   call garchx_filter_derivatives(y, spec, pars, sigma2, deriv, status, xreg, 1.5_dp)
   call assert_true(status == 0, 'derivative status')
   h = 1.0e-6_dp
   do j = 1, size(pars)
      pplus = pars
      pminus = pars
      pplus(j) = pplus(j)+h
      pminus(j) = pminus(j)-h
      call garchx_filter(y, spec, pplus, sigma_plus, status, xreg, 1.5_dp)
      call assert_true(status == 0, 'plus filter')
      call garchx_filter(y, spec, pminus, sigma_minus, status, xreg, 1.5_dp)
      call assert_true(status == 0, 'minus filter')
      fd = (sigma_plus-sigma_minus)/(2.0_dp*h)
      call assert_close(maxval(abs(fd(3:)-deriv(3:, j))), 0.0_dp, 2.0e-6_dp, &
                        'analytic derivatives')
   end do

   call garchx_objective_value(y, spec, pars, objective, status, xreg, 1.5_dp)
   objective_manual = sum(y(4:)**2/expected(4:)+log(expected(4:)))/5.0_dp
   call assert_true(status == 0, 'objective status')
   call assert_close(objective, objective_manual, 2.0e-12_dp, 'objective value')

   z = [0.2_dp, -0.4_dp, 0.6_dp, -0.8_dp, 0.1_dp, -0.2_dp, 0.3_dp, -0.5_dp]
   call garchx_simulate(8, spec, pars, ysim, ssim, zout, status, xreg, z)
   call assert_true(status == 0, 'simulation status')
   call assert_close(maxval(abs(zout-z)), 0.0_dp, 1.0e-15_dp, 'supplied innovations')
   call assert_close(maxval(abs(ysim-sqrt(ssim)*z)), 0.0_dp, 1.0e-12_dp, 'simulation identity')
   call assert_true(all(ssim > 0.0_dp), 'positive simulation variances')

   call glag(y, 2, lagged, status, pad=.false.)
   call assert_true(status == 0 .and. size(lagged) == 6, 'glag dimensions')
   call assert_close(maxval(abs(lagged-y(:6))), 0.0_dp, 1.0e-15_dp, 'glag values')
   call gdiff(y, 2, differenced, status, pad=.false.)
   call assert_true(status == 0 .and. size(differenced) == 6, 'gdiff dimensions')
   call assert_close(maxval(abs(differenced-(y(3:)-y(:6)))), 0.0_dp, 1.0e-15_dp, 'gdiff values')

   call assert_close(normal_cdf(0.0_dp), 0.5_dp, 1.0e-15_dp, 'normal cdf')
   call assert_close(normal_quantile(0.975_dp), 1.95996398454005_dp, 5.0e-9_dp, 'normal quantile')
   call assert_close(student_t_cdf(student_t_quantile(0.95_dp, 10.0_dp), 10.0_dp), &
                     0.95_dp, 2.0e-11_dp, 'student t inversion')

   matrix_data = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                          5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp], [4, 2])
   call glag(matrix_data, 1, matrix_lag, status, pad=.false.)
   call assert_true(status == 0 .and. size(matrix_lag, 1) == 3, 'matrix glag dimensions')
   call assert_close(maxval(abs(matrix_lag-matrix_data(:3, :))), 0.0_dp, 1.0e-15_dp, &
                     'matrix glag values')
   call gdiff(matrix_data, 1, matrix_diff, status, pad=.false.)
   call assert_true(status == 0 .and. size(matrix_diff, 1) == 3, 'matrix gdiff dimensions')
   call assert_close(maxval(abs(matrix_diff-(matrix_data(2:, :)-matrix_data(:3, :)))), &
                     0.0_dp, 1.0e-15_dp, 'matrix gdiff values')

   yzero = y
   yzero(5) = 0.0_dp
   call garchx_objective_value(yzero, spec, pars, objective, status, xreg, 1.5_dp, &
                               objective_mode=0)
   call assert_true(status == 0, 'zero-robust objective status')

   print '(a)', 'Core recursion, derivative, simulation, utility, and probability tests passed.'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected_value, tolerance, message)
      real(dp), intent(in) :: actual, expected_value, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual-expected_value) > tolerance) then
         write(*, '(a,2(1x,es18.10))') 'FAILED: '//trim(message), actual, expected_value
         error stop 1
      end if
   end subroutine assert_close
end program test_core
