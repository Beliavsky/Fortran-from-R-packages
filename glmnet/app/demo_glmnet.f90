! SPDX-License-Identifier: GPL-2.0-only
program demo_glmnet
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, &
      glmnet_cv_result, fit_glmnet, cv_glmnet, glmnet_family_gaussian
   implicit none
   integer, parameter :: n = 75, p = 6
   real(dp) :: x(n, p), y(n)
   integer :: i, j
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit
   type(glmnet_cv_result) :: cv
   do j = 1, p
      do i = 1, n
         x(i, j) = sin(0.03_dp * real(i * (j + 1), dp)) + &
            0.2_dp * cos(0.07_dp * real(i + 3 * j, dp))
      end do
   end do
   y = 1.0_dp + 1.5_dp * x(:, 1) - 0.9_dp * x(:, 3) + 0.5_dp * x(:, 6)
   control%alpha = 0.7_dp
   control%nlambda = 20
   call fit_glmnet(x, y, 'gaussian', fit, control)
   call cv_glmnet(x, y, glmnet_family_gaussian, cv, control, nfolds=5, seed=101)
   write(*, '(a,i0)') 'fitted lambda values: ', fit%nlambda
   write(*, '(a,es12.4)') 'cross-validated lambda.min: ', cv%lambda_min
   write(*, '(a,*(f9.4,1x))') 'coefficients at lambda.min path index: ', &
      fit%beta(:, 1, min(cv%index_min, fit%nlambda))
end program demo_glmnet
