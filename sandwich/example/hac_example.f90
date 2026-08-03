! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program hac_example
   use sandwich, only : dp, ols_model, fit_ols, newey_west_weights, vcov_hac, &
      hac_diagnostics, SANDWICH_SUCCESS
   implicit none

   integer, parameter :: n = 80
   real(dp) :: x(n, 2), y(n), noise
   real(dp), allocatable :: weights(:), covariance(:, :)
   type(ols_model) :: model
   type(hac_diagnostics) :: diagnostics
   integer :: i, status

   noise = 0.0_dp
   do i = 1, n
      x(i, :) = [1.0_dp, sin(0.08_dp * real(i, dp))]
      noise = 0.65_dp * noise + 0.2_dp * cos(0.47_dp * real(i, dp))
      y(i) = 1.0_dp + 0.5_dp * x(i, 2) + noise
   end do
   call fit_ols(x, y, model, status)
   if (status /= SANDWICH_SUCCESS) error stop
   call newey_west_weights(model%scores, weights, status, prewhite_order = 0)
   if (status /= SANDWICH_SUCCESS) error stop
   call vcov_hac(model%scores, model%bread, weights, covariance, status, &
      adjust = .false., prewhite_order = 0, diagnostics = diagnostics)
   if (status /= SANDWICH_SUCCESS) error stop

   print '(a,i0)', 'retained lag weights: ', size(weights)
   print '(a,2f12.6)', 'HAC variances: ', covariance(1, 1), covariance(2, 2)
end program hac_example
