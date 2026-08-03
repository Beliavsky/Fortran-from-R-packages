! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program hc_covariance_example
   use sandwich, only : dp, ols_model, fit_ols, vcov_hc, SANDWICH_SUCCESS
   implicit none

   real(dp) :: x(8, 2), y(8)
   real(dp), allocatable :: robust(:, :)
   type(ols_model) :: model
   integer :: i, status

   do i = 1, 8
      x(i, :) = [1.0_dp, real(i - 4, dp)]
      y(i) = 2.0_dp + 0.8_dp * x(i, 2) + 0.1_dp * real(i - 4, dp)**2 * (-1.0_dp)**i
   end do
   call fit_ols(x, y, model, status)
   if (status /= SANDWICH_SUCCESS) error stop
   call vcov_hc(x, model%residuals, model%bread, 'HC3', robust, status, model%hat)
   if (status /= SANDWICH_SUCCESS) error stop

   print '(a,2f12.6)', 'beta = ', model%coefficients
   print '(a,2f12.6)', 'HC3 variances = ', robust(1, 1), robust(2, 2)
end program hc_covariance_example
