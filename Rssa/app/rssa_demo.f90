! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program rssa_demo
   use rssa, only : dp, reconstruct_ssa, rforecast_ssa, ssa_1d, ssa_result
   implicit none
   real(dp) :: x(8)
   real(dp), allocatable :: reconstruction(:), prediction(:)
   type(ssa_result) :: fit
   integer :: i
   x = [(2.0_dp**real(i - 1, dp), i=1, 8)]
   fit = ssa_1d(x, 4, 1)
   reconstruction = reconstruct_ssa(fit, [1])
   prediction = rforecast_ssa(fit, [1], 2, only_new=.true.)
   print '(a,es14.6)', 'leading singular value: ', fit%sigma(1)
   print '(a,es14.6)', 'reconstruction max error: ', maxval(abs(reconstruction - x))
   print '(a,2es14.6)', 'next values: ', prediction
end program rssa_demo
