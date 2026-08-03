! SPDX-License-Identifier: GPL-2.0-only
program classification_example
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, &
      glmnet_assessment_result, fit_glmnet, assess_glmnet
   implicit none
   integer, parameter :: n = 120, p = 4
   real(dp) :: x(n, p), y(n), probability(n)
   integer :: i
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit
   type(glmnet_assessment_result) :: assessment
   do i = 1, n
      x(i, 1) = sin(0.09_dp * real(i, dp))
      x(i, 2) = cos(0.13_dp * real(i, dp))
      x(i, 3) = real(mod(23 * i, 41) - 20, dp) / 20.0_dp
      x(i, 4) = real(mod(11 * i, 37) - 18, dp) / 18.0_dp
   end do
   probability = 1.0_dp / (1.0_dp + exp(-(0.2_dp + 1.2_dp * x(:, 1) - 0.7_dp * x(:, 2))))
   do i = 1, n
      y(i) = merge(1.0_dp, 0.0_dp, real(mod(47 * i, 127), dp) / 127.0_dp < probability(i))
   end do
   control%alpha = 1.0_dp
   control%nlambda = 25
   control%threshold = 1.0e-6_dp
   call fit_glmnet(x, y, 'binomial', fit, control)
   call assess_glmnet(fit, x, y, assessment, 'auc')
   write(*, '(a,f8.5)') 'largest training AUC: ', maxval(assessment%value)
end program classification_example
