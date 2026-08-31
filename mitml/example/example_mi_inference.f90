! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program example_mi_inference
   use mitml, only : dp, pool_estimates, pooled_estimates
   implicit none
   real(dp) :: estimates(2, 3)
   real(dp) :: covariance(2, 2, 3)
   type(pooled_estimates) :: pooled
   integer :: i

   estimates(:, 1) = [0.40_dp, -0.20_dp]
   estimates(:, 2) = [0.44_dp, -0.16_dp]
   estimates(:, 3) = [0.39_dp, -0.19_dp]
   do i = 1, 3
      covariance(:, :, i) = reshape([0.010_dp, 0.001_dp, 0.001_dp, 0.016_dp], [2, 2])
   end do
   call pool_estimates(estimates, pooled, covariance, 200.0_dp)
   write (*, '(a,2(f10.5,1x))') "Pooled estimates: ", pooled%estimate
   write (*, '(a,2(f10.5,1x))') "Standard errors:  ", pooled%std_error
end program example_mi_inference
