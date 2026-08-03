! SPDX-License-Identifier: GPL-2.0-only
program gaussian_path_example
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, fit_glmnet
   implicit none
   integer, parameter :: n = 100, p = 5
   real(dp) :: x(n, p), y(n)
   integer :: i
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit
   do i = 1, n
      x(i, 1) = sin(0.07_dp * real(i, dp))
      x(i, 2) = cos(0.11_dp * real(i, dp))
      x(i, 3) = real(mod(17 * i, 31) - 15, dp) / 15.0_dp
      x(i, 4) = real(mod(7 * i, 19) - 9, dp) / 9.0_dp
      x(i, 5) = real(mod(13 * i, 29) - 14, dp) / 14.0_dp
   end do
   y = 0.8_dp + 1.7_dp * x(:, 1) - 1.1_dp * x(:, 2) + 0.4_dp * x(:, 4)
   control%alpha = 0.8_dp
   control%nlambda = 30
   call fit_glmnet(x, y, 'gaussian', fit, control)
   write(*, '(a,i0)') 'path length: ', fit%nlambda
   write(*, '(a,*(f10.5,1x))') 'last coefficients: ', fit%beta(:, 1, fit%nlambda)
   write(*, '(a,f10.6)') 'last deviance ratio: ', fit%dev_ratio(fit%nlambda)
end program gaussian_path_example
