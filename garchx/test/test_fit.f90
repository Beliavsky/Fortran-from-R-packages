! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
program test_fit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchx_kinds, only : dp
   use garchx_math, only : set_random_seed
   use garchx_model, only : garchx_spec, garchx_fit, make_garchx_spec
   use garchx_model, only : garchx_simulate, fit_garchx, garchx_covariance
   use garchx_model, only : garchx_forecast, garchx_quantile_path, refit_garchx
   use garchx_model, only : garchx_asymptotic_covariance, garchx_max_lag
   implicit none
   integer, parameter :: n = 900
   type(garchx_spec) :: spec
   type(garchx_fit) :: fit, new_fit, reestimated_fit
   integer :: status, i
   real(dp) :: true_par(5), initial(5), lower(5), upper(5)
   real(dp) :: xreg(n, 1), future_xreg(3, 1), new_xreg(120, 1)
   real(dp), allocatable :: avar_xreg(:, :)
   real(dp), allocatable :: y(:), sigma2(:), z(:), robust(:, :), hac(:, :)
   real(dp), allocatable :: forecast(:), paths(:, :), quantiles(:, :), avar(:, :)
   real(dp), allocatable :: new_y(:), new_sigma2(:), new_z(:)

   call set_random_seed(42017)
   call make_garchx_spec(spec, arch_lags=[1], garch_lags=[1], asym_lags=[1], xreg_count=1)
   do i = 1, n
      xreg(i, 1) = 0.15_dp + 0.05_dp*sin(0.03_dp*real(i, dp))
   end do
   true_par = [0.12_dp, 0.07_dp, 0.82_dp, 0.04_dp, 0.18_dp]
   call garchx_simulate(n, spec, true_par, y, sigma2, z, status, xreg)
   call assert_true(status == 0, 'data simulation')
   initial = [0.15_dp, 0.1_dp, 0.75_dp, 0.03_dp, 0.1_dp]
   lower = 0.0_dp
   upper = [2.0_dp, 0.8_dp, 0.99_dp, 0.8_dp, 2.0_dp]
   call fit_garchx(y, spec, fit, xreg, initial, lower, upper, vcov_type='ordinary', &
                   max_iter=2500, rel_tol=2.0e-8_dp)
   call assert_true(fit%status == 0, 'fit convergence')
   call assert_true(all(ieee_is_finite(fit%par)), 'finite parameters')
   call assert_true(all(fit%sigma2 > 0.0_dp), 'positive fitted variances')
   call assert_true(all(ieee_is_finite(fit%hessian)), 'finite hessian')
   call assert_true(all(ieee_is_finite(fit%vcov)), 'finite ordinary covariance')
   call assert_true(maxval(abs(fit%vcov-transpose(fit%vcov))) < 2.0e-8_dp, &
                    'ordinary covariance symmetry')
   call assert_true(abs(fit%par(3)-true_par(3)) < 0.20_dp, 'garch coefficient recovery')

   call garchx_covariance(y, spec, fit%par, fit%hessian, 'robust', robust, status, &
                          xreg, fit%backcast, fit%objective_mode)
   call assert_true(status == 0 .and. all(ieee_is_finite(robust)), 'robust covariance')
   call assert_true(maxval(abs(robust-transpose(robust))) < 2.0e-8_dp, &
                    'robust covariance symmetry')
   call garchx_covariance(y, spec, fit%par, fit%hessian, 'hac', hac, status, &
                          xreg, fit%backcast, fit%objective_mode, bandwidth=4.5_dp)
   call assert_true(status == 0 .and. all(ieee_is_finite(hac)), 'hac covariance')
   call assert_true(maxval(abs(hac-transpose(hac))) < 2.0e-8_dp, 'hac covariance symmetry')

   future_xreg(:, 1) = [0.14_dp, 0.15_dp, 0.16_dp]
   call garchx_forecast(fit, 3, forecast, status, future_xreg, n_sim=400, paths=paths)
   call assert_true(status == 0, 'forecast status')
   call assert_true(all(forecast > 0.0_dp) .and. all(ieee_is_finite(paths)), 'forecast values')
   call assert_true(maxval(abs(forecast-sum(paths, dim=2)/real(size(paths, 2), dp))) < 1.0e-12_dp, &
                    'forecast path average')

   call garchx_quantile_path(fit, [0.025_dp, 0.5_dp, 0.975_dp], quantiles)
   call assert_true(size(quantiles, 1) == n-garchx_max_lag(spec), 'quantile rows')
   call assert_true(all(ieee_is_finite(quantiles)), 'quantile finite')

   do i = 1, 120
      new_xreg(i, 1) = 0.13_dp + 0.03_dp*cos(0.05_dp*real(i, dp))
   end do
   call garchx_simulate(120, spec, true_par, new_y, new_sigma2, new_z, status, new_xreg)
   call assert_true(status == 0, 'new data simulation')
   call refit_garchx(fit, new_y, new_fit, status, new_xreg, reestimate=.false.)
   call assert_true(status == 0, 'fixed-parameter refit')
   call assert_true(maxval(abs(new_fit%par-fit%par)) < 1.0e-15_dp, 'refit coefficients')
   call assert_true(all(new_fit%sigma2 > 0.0_dp), 'refit variances')

   allocate(avar_xreg(2500, 1))
   do i = 1, 2500
      avar_xreg(i, 1) = 0.15_dp + 0.05_dp*sin(0.03_dp*real(i, dp))
   end do
   call garchx_asymptotic_covariance(true_par, spec, 2500, avar, status, &
                                      xreg=avar_xreg, e_eta4=3.0_dp)
   call assert_true(status == 0 .and. all(ieee_is_finite(avar)), &
                    'asymptotic covariance simulation')
   call assert_true(maxval(abs(avar-transpose(avar))) < 2.0e-8_dp, &
                    'asymptotic covariance symmetry')

   call refit_garchx(fit, y, reestimated_fit, status, xreg, reestimate=.true., &
                     max_iter=800, rel_tol=1.0e-7_dp)
   call assert_true(status == 0, 'reestimated refit')
   call assert_true(all(ieee_is_finite(reestimated_fit%par)), 'reestimated parameters')

   print '(a)', 'QML fitting, covariance, forecasting, quantile, refit, and avar tests passed.'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_fit
