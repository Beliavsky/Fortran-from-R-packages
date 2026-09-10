! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program forecast_gapfill
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use rssa, only : dp, igapfill, rforecast_ssa, ssa_1d, ssa_result
   implicit none
   real(dp) :: x(8), incomplete(8)
   real(dp), allocatable :: filled(:), prediction(:)
   type(ssa_result) :: fit
   integer :: i, info
   x = [(2.0_dp**real(i - 1, dp), i=1, 8)]
   incomplete = x
   incomplete(4) = ieee_value(0.0_dp, ieee_quiet_nan)
   filled = igapfill(incomplete, 4, 1, maxiter=50, info=info)
   fit = ssa_1d(filled, 4, 1)
   prediction = rforecast_ssa(fit, [1], 2, only_new=.true., info=info)
   print '(a,es14.6)', 'filled fourth value: ', filled(4)
   print '(a,2es14.6)', 'forecast: ', prediction
end program forecast_gapfill
