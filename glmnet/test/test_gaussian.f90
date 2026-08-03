! SPDX-License-Identifier: GPL-2.0-only
program test_gaussian
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, glmnet_cv_result, &
      glmnet_sparse_csc, fit_glmnet, fit_glmnet_sparse, predict_glmnet, &
      predict_glmnet_at, coef_glmnet, nonzero_coef, cv_glmnet, relax_glmnet, &
      dense_to_sparse, glmnet_family_gaussian, glmnet_success
   implicit none
   integer, parameter :: n = 80, p = 4
   real(dp) :: x(n, p), y(n), lambda(4), lower(p), upper(p)
   real(dp), allocatable :: prediction(:,:,:), prediction_at(:,:,:), a(:), b(:,:)
   integer, allocatable :: nz(:)
   integer :: i, status
   logical :: excluded(p)
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit, sparse_fit, bounded, relaxed
   type(glmnet_cv_result) :: cv
   type(glmnet_sparse_csc) :: sparse

   do i = 1, n
      x(i, 1) = sin(0.13_dp * real(i, dp))
      x(i, 2) = cos(0.07_dp * real(i, dp))
      x(i, 3) = real(mod(17 * i, 31) - 15, dp) / 15.0_dp
      x(i, 4) = real(mod(11 * i, 23) - 11, dp) / 11.0_dp
   end do
   y = 1.25_dp + 2.0_dp * x(:, 1) - 1.5_dp * x(:, 2) + 0.35_dp * x(:, 3)
   lambda = [0.5_dp, 0.1_dp, 0.01_dp, 1.0e-8_dp]
   control%alpha = 0.8_dp
   control%threshold = 1.0e-10_dp
   control%max_iterations = 20000

   call fit_glmnet(x, y, 'gaussian', fit, control, lambda=lambda)
   call check(fit%status == glmnet_success, 'gaussian status')
   call check(fit%nlambda == size(lambda), 'gaussian path length')
   call check(maxval(abs(fit%beta(:, 1, 4) - [2.0_dp, -1.5_dp, 0.35_dp, 0.0_dp])) < 2.0e-4_dp, &
      'gaussian coefficients')
   call check(abs(fit%intercept(1, 4) - 1.25_dp) < 2.0e-4_dp, 'gaussian intercept')
   call check(fit%dev_ratio(4) > 0.999999_dp, 'gaussian deviance')

   call predict_glmnet(fit, x, prediction, status)
   call check(status == glmnet_success, 'prediction status')
   call check(sqrt(sum((prediction(:, 1, 4) - y) ** 2) / real(n, dp)) < 2.0e-4_dp, &
      'prediction accuracy')
   call predict_glmnet_at(fit, x, [0.05_dp], prediction_at, status)
   call check(status == glmnet_success .and. size(prediction_at, 3) == 1, 'interpolated prediction')
   call coef_glmnet(fit, 0.1_dp, a, b, status)
   call check(status == glmnet_success .and. size(b, 1) == p, 'coefficient extraction')
   call nonzero_coef(fit, 4, nz, status)
   call check(status == glmnet_success .and. size(nz) == 3, 'nonzero coefficients')

   lower = -huge(1.0_dp)
   upper = huge(1.0_dp)
   lower(2) = 0.0_dp
   excluded = .false.
   excluded(3) = .true.
   call fit_glmnet(x, y, 'gaussian', bounded, control, lambda=lambda, lower=lower, &
      upper=upper, excluded=excluded)
   call check(abs(bounded%beta(2, 1, 4)) < 1.0e-12_dp, 'lower bound')
   call check(abs(bounded%beta(3, 1, 4)) < 1.0e-12_dp, 'excluded predictor')

   call dense_to_sparse(x, sparse, 0.0_dp)
   call fit_glmnet_sparse(sparse, y, 'gaussian', sparse_fit, control, lambda=lambda)
   call check(maxval(abs(sparse_fit%beta - fit%beta)) < 1.0e-10_dp, 'sparse wrapper')

   control%nlambda = 8
   control%lambda_min_ratio = 1.0e-3_dp
   call cv_glmnet(x, y, glmnet_family_gaussian, cv, control, nfolds=5, seed=19)
   call check(cv%status == glmnet_success, 'gaussian cv status')
   call check(cv%lambda_min > 0.0_dp .and. cv%index_1se <= cv%index_min, 'gaussian cv selection')

   call relax_glmnet(fit, x, y, relaxed, control)
   call check(relaxed%dev_ratio(2) >= fit%dev_ratio(2) - 1.0e-8_dp, 'relaxed fit')
   print '(a)', 'test_gaussian: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine check
end program test_gaussian
