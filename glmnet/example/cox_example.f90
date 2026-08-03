! SPDX-License-Identifier: GPL-2.0-only
program cox_example
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, &
      fit_cox_path, predict_glmnet, concordance_index
   implicit none
   integer, parameter :: n = 80, p = 3
   real(dp) :: x(n, p), start_time(n), stop_time(n)
   real(dp), allocatable :: prediction(:,:,:)
   integer :: event(n), i, status
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit
   do i = 1, n
      x(i, 1) = sin(0.12_dp * real(i, dp))
      x(i, 2) = cos(0.08_dp * real(i, dp))
      x(i, 3) = real(mod(19 * i, 31) - 15, dp) / 15.0_dp
      start_time(i) = 0.0_dp
      stop_time(i) = exp(-0.8_dp * x(i, 1) + 0.25_dp * x(i, 2)) * &
         (0.8_dp + real(mod(13 * i, 23), dp) / 12.0_dp)
      event(i) = merge(0, 1, mod(i, 5) == 0)
   end do
   control%alpha = 0.9_dp
   control%nlambda = 20
   control%threshold = 1.0e-6_dp
   call fit_cox_path(x, start_time, stop_time, event, fit, control)
   call predict_glmnet(fit, x, prediction, status, prediction_type='link')
   write(*, '(a,f8.5)') 'training C-index: ', &
      concordance_index(stop_time, event, prediction(:, 1, fit%nlambda))
end program cox_example
