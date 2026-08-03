! SPDX-License-Identifier: GPL-2.0-only
program test_cox
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, glmnet_cv_result, &
      fit_cox_path, predict_glmnet, concordance_index, cv_cox, glmnet_success, &
      glmnet_max_iterations
   use glmnet_cox, only : cox_loss_gradient
   implicit none
   integer, parameter :: n = 90, p = 3
   real(dp) :: x(n, p), start_time(n), stop_time(n), offset(n), lambda(3)
   real(dp) :: weights(n), beta(p), loss, loss_plus, loss_minus, h
   real(dp), allocatable :: gradient(:), dummy(:), prediction(:,:,:)
   integer :: event(n), strata(n), i, j, status
   real(dp) :: numerical, cvalue
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit, tied_fit
   type(glmnet_cv_result) :: cv

   do i = 1, n
      x(i, 1) = sin(0.12_dp * real(i, dp))
      x(i, 2) = cos(0.08_dp * real(i, dp))
      x(i, 3) = real(mod(23 * i, 41) - 20, dp) / 20.0_dp
      start_time(i) = 0.0_dp
      stop_time(i) = exp(-0.9_dp * x(i, 1) + 0.35_dp * x(i, 2)) * &
         (0.8_dp + real(mod(17 * i, 29), dp) / 15.0_dp)
      event(i) = merge(0, 1, mod(i, 5) == 0)
      strata(i) = merge(1, 2, i <= n / 2)
      weights(i) = 1.0_dp
      offset(i) = 0.0_dp
   end do
   lambda = [0.15_dp, 0.04_dp, 0.008_dp]
   control%alpha = 0.9_dp
   control%threshold = 2.0e-6_dp
   control%max_iterations = 10000

   call fit_cox_path(x, start_time, stop_time, event, fit, control, weights_in=weights, &
      offset_in=offset, strata_in=strata, lambda_in=lambda)
   call check(fit%status == glmnet_success .or. fit%status == glmnet_max_iterations, 'cox status')
   call predict_glmnet(fit, x, prediction, status, prediction_type='link')
   call check(status == glmnet_success, 'cox prediction')
   cvalue = concordance_index(stop_time, event, prediction(:, 1, 3))
   call check(cvalue > 0.62_dp, 'cox concordance')

   allocate(gradient(p), dummy(p))
   beta = [0.2_dp, -0.1_dp, 0.05_dp]
   call cox_loss_gradient(x, start_time, stop_time, event, strata, weights / sum(weights), &
      offset, beta, .false., loss, gradient, status)
   call check(status == glmnet_success, 'cox gradient status')
   h = 1.0e-6_dp
   do j = 1, p
      beta(j) = beta(j) + h
      call cox_loss_gradient(x, start_time, stop_time, event, strata, weights / sum(weights), &
         offset, beta, .false., loss_plus, dummy, status)
      beta(j) = beta(j) - 2.0_dp * h
      call cox_loss_gradient(x, start_time, stop_time, event, strata, weights / sum(weights), &
         offset, beta, .false., loss_minus, dummy, status)
      beta(j) = beta(j) + h
      numerical = (loss_plus - loss_minus) / (2.0_dp * h)
      call check(abs(numerical - gradient(j)) < 2.0e-5_dp, 'cox finite-difference gradient')
   end do

   stop_time = 1.0_dp + real([(mod(i, 9), i=1,n)], dp)
   call fit_cox_path(x, start_time, stop_time, event, tied_fit, control, lambda_in=lambda, efron=.true.)
   call check(tied_fit%nlambda == 3, 'efron tied fit')

   control%nlambda = 3
   call cv_cox(x, start_time, stop_time, event, cv, control, nfolds=3, seed=11, efron=.true.)
   call check(cv%status == glmnet_success, 'cox cv status')
   call check(cv%lambda_min > 0.0_dp, 'cox cv lambda')
   print '(a)', 'test_cox: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine check
end program test_cox
