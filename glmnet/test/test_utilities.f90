! SPDX-License-Identifier: GPL-2.0-only
program test_utilities
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, glmnet_roc_result, &
      glmnet_survival_data, fit_custom_family_path, gaussian_identity_working, &
      na_replace, prepare_x, make_x, rmult, stratify_surv, auc, roc_glmnet, &
      confusion_glmnet, glmnet_measures, glmnet_success, glmnet_invalid_argument, &
      fit_cox_path
   implicit none
   integer, parameter :: n = 40, p = 3
   real(dp) :: x(n, p), y(n), xnan(4, 3), table(2, 2), lambda(2)
   real(dp), allocatable :: filled(:,:), prepared(:,:), combined(:,:)
   integer, allocatable :: keep(:), counts(:)
   character(len=16), allocatable :: measures(:)
   integer :: i, status
   integer :: event(n), strata(n)
   real(dp) :: start_time(n), stop_time(n)
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit, invalid_fit
   type(glmnet_roc_result) :: roc
   type(glmnet_survival_data) :: survival

   do i = 1, n
      x(i, 1) = real(i, dp) / real(n, dp)
      x(i, 2) = sin(0.2_dp * real(i, dp))
      x(i, 3) = 1.0_dp
   end do
   y = 0.7_dp + 1.8_dp * x(:, 1) - 0.4_dp * x(:, 2)
   lambda = [0.1_dp, 1.0e-7_dp]
   control%threshold = 1.0e-9_dp
   call fit_custom_family_path(x, y, gaussian_identity_working, fit, control, lambda_in=lambda)
   call check(fit%status == glmnet_success, 'custom family status')
   call check(maxval(abs(fit%beta(:2, 1, 2) - [1.8_dp, -0.4_dp])) < 1.0e-3_dp, &
      'custom family coefficients')

   xnan = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
      2.0_dp, ieee_value(1.0_dp, ieee_quiet_nan), 4.0_dp, 6.0_dp, &
      5.0_dp, 5.0_dp, 5.0_dp, 5.0_dp], shape(xnan))
   call na_replace(xnan, filled, status)
   call check(status == glmnet_success .and. abs(filled(2, 2) - 4.0_dp) < 1.0e-12_dp, &
      'na replacement')
   call prepare_x(xnan, prepared, keep, status)
   call check(status == glmnet_success .and. size(prepared, 2) == 2, 'prepare x')
   call make_x(prepared, combined, status, test=prepared(:2, :))
   call check(status == glmnet_success .and. size(combined, 1) == 6, 'make x')

   call rmult([0.2_dp, 0.3_dp, 0.5_dp], 1000, counts, status, seed=3)
   call check(status == glmnet_success .and. sum(counts) == 1000, 'rmult total')
   call check(all(counts > 100), 'rmult distribution')

   y = merge(1.0_dp, 0.0_dp, x(:, 1) > 0.5_dp)
   call roc_glmnet(y, x(:, 1), roc)
   call check(roc%auc > 0.99_dp .and. abs(auc(y, x(:, 1)) - roc%auc) < 1.0e-12_dp, &
      'roc and auc')
   call confusion_glmnet(y, x(:, 1), 0.5_dp, table, status)
   call check(status == glmnet_success .and. abs(sum(table) - real(n, dp)) < 1.0e-12_dp, 'confusion table')
   call glmnet_measures(2, measures)
   call check(size(measures) == 4, 'measure names')

   start_time = 0.0_dp
   stop_time = real([(i, i=1,n)], dp)
   event = 1
   strata = 1 + mod([(i, i=1,n)], 2)
   call stratify_surv(start_time, stop_time, event, strata, survival, status)
   call check(status == glmnet_success .and. size(survival%stop) == n, 'stratified survival')
   event = 0
   call fit_cox_path(x(:, :2), start_time, stop_time, event, invalid_fit, control)
   call check(invalid_fit%status /= glmnet_success, 'no-event error')
   call stratify_surv(start_time(1:3), stop_time(1:2), event(1:3), strata(1:3), survival, status)
   call check(status == glmnet_invalid_argument, 'invalid survival dimensions')
   print '(a)', 'test_utilities: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine check
end program test_utilities
