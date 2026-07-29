! SPDX-License-Identifier: GPL-3.0-only
program mcrp_demo
   use mcrp_module
   implicit none

   real(dp) :: r(12, 3), start(3), lower(3)
   type(mcrp_result) :: fit
   integer :: i

   do i = 1, size(r, 1)
      r(i, 1) = 0.010_dp * sin(0.7_dp * real(i, dp))
      r(i, 2) = 0.007_dp * cos(0.5_dp * real(i, dp)) + 0.2_dp * r(i, 1)
      r(i, 3) = 0.005_dp * sin(0.9_dp * real(i, dp)) + &
         0.0002_dp * real(i * i, dp)
   end do
   start = 1.0_dp / 3.0_dp
   lower = 0.0_dp

   call mcrp(start, r, fit, lower=lower)
   write(*, '(a, l1)') 'Converged: ', fit%converged
   write(*, '(a, es14.6)') 'Objective: ', fit%objective
   write(*, '(a, *(f10.6, 1x))') 'Weights: ', fit%weights
   write(*, '(a, *(f10.6, 1x))') 'Variance contributions: ', &
      fit%variance_contributions
   write(*, '(a, *(f10.6, 1x))') 'Skewness contributions: ', &
      fit%skewness_contributions
   write(*, '(a, *(f10.6, 1x))') 'Kurtosis contributions: ', &
      fit%kurtosis_contributions
end program mcrp_demo
