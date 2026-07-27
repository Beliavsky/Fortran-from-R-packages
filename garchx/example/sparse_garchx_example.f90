! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
program sparse_garchx_example
   use garchx_kinds, only : dp
   use garchx_model, only : garchx_spec, make_garchx_spec, garchx_simulate
   implicit none
   integer, parameter :: n = 12
   integer :: status, i
   type(garchx_spec) :: spec
   real(dp) :: pars(6), innovations(n)
   real(dp), allocatable :: y(:), sigma2(:), z(:)

   call make_garchx_spec(spec, arch_lags=[1, 4], garch_lags=[1, 3], asym_lags=[2])
   pars = [0.15_dp, 0.08_dp, 0.02_dp, 0.60_dp, 0.10_dp, 0.04_dp]
   innovations = [(0.08_dp*real(i-6, dp), i=1, n)]
   call garchx_simulate(n, spec, pars, y, sigma2, z, status, &
                        supplied_innovations=innovations)
   if (status /= 0) error stop 'sparse simulation failed'
   print '(a)', ' t           y       sigma2'
   do i = 1, n
      print '(i2,2(1x,f12.6))', i, y(i), sigma2(i)
   end do
end program sparse_garchx_example
