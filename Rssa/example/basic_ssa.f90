! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program basic_ssa
   use rssa, only : contributions, dp, reconstruct_ssa, ssa_1d, ssa_result, wcor_ssa
   implicit none
   real(dp) :: x(12)
   real(dp), allocatable :: share(:), reconstruction(:), correlations(:, :)
   type(ssa_result) :: fit
   integer :: i
   x = [(sin(0.4_dp * real(i, dp)) + 0.1_dp * real(i, dp), i=1, 12)]
   fit = ssa_1d(x, 6, 4)
   share = contributions(fit)
   reconstruction = reconstruct_ssa(fit, [1, 2])
   correlations = wcor_ssa(fit, [1, 2, 3, 4])
   print '(a,4f10.5)', 'component contributions: ', share
   print '(a,es14.6)', 'two-component reconstruction norm: ', sqrt(sum(reconstruction**2))
   print '(a,es14.6)', 'w-correlation(1,2): ', correlations(1, 2)
end program basic_ssa
