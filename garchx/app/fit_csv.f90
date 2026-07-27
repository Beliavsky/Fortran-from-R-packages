! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
program fit_csv
   use garchx_kinds, only : dp
   use garchx_csv, only : read_dated_numeric_csv, parse_integer_list
   use garchx_model, only : garchx_spec, garchx_fit, make_garchx_spec, fit_garchx
   implicit none
   character(len=512) :: filename, arch_text, garch_text, asym_text, covariance_type
   character(len=64), allocatable :: dates(:)
   integer, allocatable :: arch_lags(:), garch_lags(:), asym_lags(:)
   real(dp), allocatable :: data(:, :), xreg(:, :)
   type(garchx_spec) :: spec
   type(garchx_fit) :: fit
   integer :: status, narg, i, kx

   narg = command_argument_count()
   if (narg < 1) then
      print '(a)', 'usage: fit_csv FILE [ARCH_LAGS] [GARCH_LAGS] [ASYM_LAGS] [VCOV]'
      print '(a)', 'example: fit_csv data/example.csv 1 1 1 hac'
      stop 2
   end if
   call get_command_argument(1, filename)
   arch_text = '1'
   garch_text = '1'
   asym_text = '-'
   covariance_type = 'ordinary'
   if (narg >= 2) call get_command_argument(2, arch_text)
   if (narg >= 3) call get_command_argument(3, garch_text)
   if (narg >= 4) call get_command_argument(4, asym_text)
   if (narg >= 5) call get_command_argument(5, covariance_type)
   call parse_integer_list(arch_text, arch_lags, status)
   if (status /= 0) error stop 'invalid ARCH lag list'
   call parse_integer_list(garch_text, garch_lags, status)
   if (status /= 0) error stop 'invalid GARCH lag list'
   call parse_integer_list(asym_text, asym_lags, status)
   if (status /= 0) error stop 'invalid ASYM lag list'
   call read_dated_numeric_csv(trim(filename), dates, data, status)
   if (status /= 0) error stop 'could not read CSV'
   kx = size(data, 2)-1
   call make_garchx_spec(spec, arch_lags, garch_lags, asym_lags, kx)
   if (kx > 0) then
      allocate(xreg(size(data, 1), kx))
      xreg = data(:, 2:)
      call fit_garchx(data(:, 1), spec, fit, xreg=xreg, vcov_type=trim(covariance_type), &
                      max_iter=3000, rel_tol=2.0e-8_dp)
   else
      call fit_garchx(data(:, 1), spec, fit, vcov_type=trim(covariance_type), &
                      max_iter=3000, rel_tol=2.0e-8_dp)
   end if
   if (fit%status /= 0) error stop 'model fitting failed'
   print '(a,1x,a,1x,a)', 'Sample:', trim(dates(1)), trim(dates(size(dates)))
   print '(a,i0)', 'Observations: ', size(data, 1)
   print '(a,1x,a)', 'Covariance type:', trim(fit%vcov_type)
   do i = 1, size(fit%par)
      print '(a,i0,a,f14.7,a,f14.7)', 'parameter(', i, ') = ', fit%par(i), &
            '  std.error = ', sqrt(max(0.0_dp, fit%vcov(i, i)))
   end do
   print '(a,f16.6)', 'Gaussian log likelihood: ', fit%loglik
end program fit_csv
