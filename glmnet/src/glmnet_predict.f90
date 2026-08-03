! SPDX-License-Identifier: GPL-2.0-only
module glmnet_predict
   use glmnet_kinds, only : dp, glmnet_eps
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument
   use glmnet_types, only : glmnet_path_result, glmnet_family_gaussian, &
      glmnet_family_binomial, glmnet_family_poisson, glmnet_family_multinomial, &
      glmnet_family_mgaussian, glmnet_family_cox
   use glmnet_utils, only : linear_interpolate_coefficients, logistic, safe_exp
   implicit none
   private
   public :: predict_glmnet, predict_glmnet_at, coef_glmnet, nonzero_coef
contains
   subroutine predict_glmnet(fit, x, prediction, status, prediction_type, offset)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: prediction(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(in), optional :: prediction_type
      real(dp), intent(in), optional :: offset(:,:)
      character(len=16) :: ptype
      real(dp), allocatable :: eta(:,:)
      integer :: n, l
      status = glmnet_success
      n = size(x, 1)
      if (size(x, 2) /= fit%nvars .or. fit%nlambda < 1) then
         status = glmnet_invalid_argument
         allocate(prediction(0, 0, 0))
         return
      end if
      ptype = 'response'
      if (present(prediction_type)) ptype = adjustl(prediction_type)
      allocate(prediction(n, fit%nout, fit%nlambda), eta(n, fit%nout))
      do l = 1, fit%nlambda
         eta = matmul(x, fit%beta(:, :, l)) + spread(fit%intercept(:, l), 1, n)
         if (present(offset)) then
            if (size(offset, 1) /= n .or. size(offset, 2) /= fit%nout) then
               status = glmnet_invalid_argument
               prediction = 0.0_dp
               return
            end if
            eta = eta + offset
         end if
         call transform_prediction(fit%family_code, eta, ptype, prediction(:, :, l))
      end do
   end subroutine predict_glmnet

   subroutine predict_glmnet_at(fit, x, s, prediction, status, prediction_type, offset)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: x(:,:), s(:)
      real(dp), allocatable, intent(out) :: prediction(:,:,:)
      integer, intent(out) :: status
      character(len=*), intent(in), optional :: prediction_type
      real(dp), intent(in), optional :: offset(:,:)
      character(len=16) :: ptype
      real(dp), allocatable :: a(:), b(:,:), eta(:,:)
      integer :: n, m
      status = glmnet_success
      n = size(x, 1)
      if (size(x, 2) /= fit%nvars .or. fit%nlambda < 1 .or. size(s) < 1) then
         status = glmnet_invalid_argument
         allocate(prediction(0, 0, 0))
         return
      end if
      ptype = 'response'
      if (present(prediction_type)) ptype = adjustl(prediction_type)
      allocate(prediction(n, fit%nout, size(s)), a(fit%nout), &
         b(fit%nvars, fit%nout), eta(n, fit%nout))
      do m = 1, size(s)
         call linear_interpolate_coefficients(fit%lambda, fit%intercept, fit%beta, &
            s(m), a, b)
         eta = matmul(x, b) + spread(a, 1, n)
         if (present(offset)) then
            if (size(offset, 1) /= n .or. size(offset, 2) /= fit%nout) then
               status = glmnet_invalid_argument
               prediction = 0.0_dp
               return
            end if
            eta = eta + offset
         end if
         call transform_prediction(fit%family_code, eta, ptype, prediction(:, :, m))
      end do
   end subroutine predict_glmnet_at

   subroutine coef_glmnet(fit, s, intercept, beta, status)
      type(glmnet_path_result), intent(in) :: fit
      real(dp), intent(in) :: s
      real(dp), allocatable, intent(out) :: intercept(:), beta(:,:)
      integer, intent(out) :: status
      status = glmnet_success
      if (fit%nlambda < 1) then
         status = glmnet_invalid_argument
         allocate(intercept(0), beta(0, 0))
         return
      end if
      allocate(intercept(fit%nout), beta(fit%nvars, fit%nout))
      call linear_interpolate_coefficients(fit%lambda, fit%intercept, fit%beta, &
         s, intercept, beta)
   end subroutine coef_glmnet

   subroutine nonzero_coef(fit, lambda_index, indices, status, tolerance)
      type(glmnet_path_result), intent(in) :: fit
      integer, intent(in) :: lambda_index
      integer, allocatable, intent(out) :: indices(:)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tolerance
      real(dp) :: tol
      integer :: j, count_nonzero
      status = glmnet_success
      tol = 1.0e-12_dp
      if (present(tolerance)) tol = max(tolerance, 0.0_dp)
      if (lambda_index < 1 .or. lambda_index > fit%nlambda) then
         status = glmnet_invalid_argument
         allocate(indices(0))
         return
      end if
      count_nonzero = 0
      do j = 1, fit%nvars
         if (maxval(abs(fit%beta(j, :, lambda_index))) > tol) count_nonzero = count_nonzero + 1
      end do
      allocate(indices(count_nonzero))
      count_nonzero = 0
      do j = 1, fit%nvars
         if (maxval(abs(fit%beta(j, :, lambda_index))) > tol) then
            count_nonzero = count_nonzero + 1
            indices(count_nonzero) = j
         end if
      end do
   end subroutine nonzero_coef

   subroutine transform_prediction(family_code, eta, ptype, prediction)
      integer, intent(in) :: family_code
      real(dp), intent(in) :: eta(:,:)
      character(len=*), intent(in) :: ptype
      real(dp), intent(out) :: prediction(size(eta, 1), size(eta, 2))
      real(dp) :: maximum, total
      integer :: i, best
      if (trim(ptype) == 'link') then
         prediction = eta
         return
      end if
      select case (family_code)
      case (glmnet_family_gaussian, glmnet_family_mgaussian)
         prediction = eta
      case (glmnet_family_binomial)
         prediction = logistic(eta)
         if (trim(ptype) == 'class') then
            where (prediction >= 0.5_dp)
               prediction = 1.0_dp
            elsewhere
               prediction = 0.0_dp
            end where
         end if
      case (glmnet_family_poisson)
         prediction = safe_exp(eta)
      case (glmnet_family_multinomial)
         do i = 1, size(eta, 1)
            maximum = maxval(eta(i, :))
            prediction(i, :) = exp(max(eta(i, :) - maximum, -700.0_dp))
            total = sum(prediction(i, :))
            prediction(i, :) = prediction(i, :) / max(total, glmnet_eps)
            if (trim(ptype) == 'class') then
               best = maxloc(prediction(i, :), dim=1)
               prediction(i, :) = 0.0_dp
               prediction(i, best) = real(best, dp)
            end if
         end do
      case (glmnet_family_cox)
         if (trim(ptype) == 'response' .or. trim(ptype) == 'risk') then
            prediction = safe_exp(eta)
         else
            prediction = eta
         end if
      case default
         prediction = eta
      end select
   end subroutine transform_prediction
end module glmnet_predict
