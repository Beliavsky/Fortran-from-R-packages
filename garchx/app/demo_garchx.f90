! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
program demo_garchx
   use garchx_kinds, only : dp
   use garchx_math, only : set_random_seed
   use garchx_model, only : garchx_spec, garchx_fit, make_garchx_spec
   use garchx_model, only : garchx_simulate, fit_garchx, garchx_forecast
   implicit none
   integer, parameter :: n = 1000
   integer :: i, status
   type(garchx_spec) :: spec
   type(garchx_fit) :: fit
   real(dp) :: true_par(5), initial(5), lower(5), upper(5), xreg(n, 1), future_xreg(5, 1)
   real(dp), allocatable :: y(:), sigma2(:), z(:), forecast(:)

   call set_random_seed(12345)
   call make_garchx_spec(spec, arch_lags=[1], garch_lags=[1], asym_lags=[1], xreg_count=1)
   do i = 1, n
      xreg(i, 1) = 0.12_dp+0.04_dp*sin(0.02_dp*real(i, dp))
   end do
   true_par = [0.12_dp, 0.08_dp, 0.80_dp, 0.05_dp, 0.20_dp]
   call garchx_simulate(n, spec, true_par, y, sigma2, z, status, xreg)
   if (status /= 0) error stop 'simulation failed'

   initial = [0.15_dp, 0.10_dp, 0.72_dp, 0.03_dp, 0.10_dp]
   lower = 0.0_dp
   upper = [2.0_dp, 0.9_dp, 0.99_dp, 0.9_dp, 2.0_dp]
   call fit_garchx(y, spec, fit, xreg, initial, lower, upper, vcov_type='hac', &
                   bandwidth=4.0_dp, max_iter=2500, rel_tol=2.0e-8_dp)
   if (fit%status /= 0) error stop 'fit failed'
   future_xreg(:, 1) = [0.11_dp, 0.12_dp, 0.13_dp, 0.14_dp, 0.15_dp]
   call garchx_forecast(fit, 5, forecast, status, future_xreg, n_sim=1000)
   if (status /= 0) error stop 'forecast failed'

   print '(a)', 'GARCH-X demonstration'
   print '(a,*(1x,f10.5))', 'True parameters:     ', true_par
   print '(a,*(1x,f10.5))', 'Estimated parameters:', fit%par
   print '(a,f14.5)', 'Gaussian log likelihood: ', fit%loglik
   print '(a,*(1x,f10.5))', 'Variance forecasts:  ', forecast
end program demo_garchx
