! SPDX-License-Identifier: GPL-2.0-or-later
program apca_example
   use fints
   implicit none

   integer, parameter :: observations = 36, series = 40
   real(dp) :: x(observations, series), factor, loading
   type(apca_result) :: fit
   integer :: i, j

   do j = 1, series
      loading = 0.20_dp + 0.035_dp * real(j, dp)
      do i = 1, observations
         factor = sin(0.18_dp * real(i, dp)) + 0.3_dp * cos(0.41_dp * real(i, dp))
         x(i, j) = loading * factor + 0.04_dp * sin(0.13_dp * real(i * j, dp))
      end do
   end do

   call apca(x, 1, fit)
   print '(a,f12.5)', 'leading eigenvalue = ', fit%eigenvalues(1)
   print '(a,f12.5)', 'mean R-squared = ', sum(fit%r_squared) / real(series, dp)
   print '(a,5f10.5)', 'first five loadings = ', fit%loadings(1:5, 1)
end program apca_example
