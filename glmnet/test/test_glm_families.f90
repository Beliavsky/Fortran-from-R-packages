! SPDX-License-Identifier: GPL-2.0-only
program test_glm_families
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, &
      glmnet_assessment_result, glmnet_cv_result, fit_glmnet, fit_multinomial_path, &
      fit_multinomial_matrix_path, predict_glmnet, assess_glmnet, assess_multinomial, &
      cv_multinomial, glmnet_success, glmnet_max_iterations
   implicit none
   integer, parameter :: n = 120, p = 4, k = 3
   real(dp) :: x(n, p), yb(n), yp(n), ym(n, 2), probabilities(n)
   real(dp) :: ymat(n, k), lambda(3)
   real(dp), allocatable :: prediction(:,:,:)
   integer :: class_id(n), i, status, correct
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: binfit, poifit, mgfit, multifit, multimatrix
   type(glmnet_assessment_result) :: assessment
   type(glmnet_cv_result) :: cv

   do i = 1, n
      x(i, 1) = sin(0.09_dp * real(i, dp))
      x(i, 2) = cos(0.14_dp * real(i, dp))
      x(i, 3) = real(mod(29 * i, 47) - 23, dp) / 23.0_dp
      x(i, 4) = real(mod(13 * i, 37) - 18, dp) / 18.0_dp
   end do
   probabilities = 1.0_dp / (1.0_dp + exp(-(0.15_dp + 1.1_dp * x(:, 1) - 0.8_dp * x(:, 2))))
   do i = 1, n
      yb(i) = merge(1.0_dp, 0.0_dp, real(mod(37 * i, 101), dp) / 101.0_dp < probabilities(i))
      yp(i) = real(max(0, nint(exp(0.25_dp + 0.55_dp * x(i, 1) - 0.25_dp * x(i, 3)) + &
         0.35_dp * sin(real(i, dp)))), dp)
   end do
   ym(:, 1) = 0.5_dp + 1.4_dp * x(:, 1) - 0.6_dp * x(:, 2)
   ym(:, 2) = -0.2_dp + 0.3_dp * x(:, 1) + 1.1_dp * x(:, 3)
   do i = 1, n
      if (x(i, 1) + 0.25_dp * x(i, 3) > 0.45_dp) then
         class_id(i) = 1
      else if (x(i, 2) - 0.2_dp * x(i, 3) > 0.25_dp) then
         class_id(i) = 2
      else
         class_id(i) = 3
      end if
      if (mod(i, 19) == 0) class_id(i) = 1 + mod(class_id(i), k)
   end do
   ymat = 0.0_dp
   do i = 1, n
      ymat(i, class_id(i)) = 1.0_dp
   end do

   lambda = [0.2_dp, 0.04_dp, 0.005_dp]
   control%alpha = 0.8_dp
   control%threshold = 2.0e-6_dp
   control%max_iterations = 10000
   control%max_outer_iterations = 100

   call fit_glmnet(x, yb, 'binomial', binfit, control, lambda=lambda)
   call check(binfit%status == glmnet_success, 'binomial status')
   call check(binfit%beta(1, 1, 3) > 0.0_dp .and. binfit%beta(2, 1, 3) < 0.0_dp, &
      'binomial signs')
   call assess_glmnet(binfit, x, yb, assessment, 'auc')
   call check(assessment%value(3) > 0.70_dp, 'binomial auc')

   call fit_glmnet(x, yp, 'poisson', poifit, control, lambda=lambda)
   call check(poifit%status == glmnet_success, 'poisson status')
   call check(poifit%beta(1, 1, 3) > 0.0_dp, 'poisson sign')
   call check(poifit%dev_ratio(3) > 0.10_dp, 'poisson deviance')

   call fit_glmnet(x, ym, 'mgaussian', mgfit, control, lambda=lambda)
   call check(mgfit%status == glmnet_success, 'mgaussian status')
   call check(maxval(abs(mgfit%beta(:, 1, 3) - [1.4_dp, -0.6_dp, 0.0_dp, 0.0_dp])) < 0.08_dp, &
      'mgaussian response 1')
   call check(maxval(abs(mgfit%beta(:, 2, 3) - [0.3_dp, 0.0_dp, 1.1_dp, 0.0_dp])) < 0.08_dp, &
      'mgaussian response 2')

   control%grouped = .true.
   control%threshold = 1.0e-5_dp
   call fit_multinomial_path(x, class_id, multifit, control, lambda_in=lambda, nclass=k)
   call check(multifit%status == glmnet_success .or. multifit%status == glmnet_max_iterations, &
      'multinomial status')
   call predict_glmnet(multifit, x, prediction, status)
   call check(status == glmnet_success, 'multinomial prediction')
   correct = 0
   do i = 1, n
      if (maxloc(prediction(i, :, 3), dim=1) == class_id(i)) correct = correct + 1
   end do
   call check(real(correct, dp) / real(n, dp) > 0.78_dp, 'multinomial accuracy')
   call fit_multinomial_matrix_path(x, ymat, multimatrix, control, lambda_in=lambda)
   call check(maxval(abs(multimatrix%beta - multifit%beta)) < 1.0e-8_dp, &
      'multinomial matrix interface')
   call assess_multinomial(multifit, x, class_id, assessment, 'class')
   call check(assessment%value(3) < 0.22_dp, 'multinomial assessment')

   control%nlambda = 3
   call cv_multinomial(x, class_id, cv, control, nfolds=3, seed=7, nclass=k)
   call check(cv%status == glmnet_success, 'multinomial cv')
   print '(a)', 'test_glm_families: PASS'
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine check
end program test_glm_families
