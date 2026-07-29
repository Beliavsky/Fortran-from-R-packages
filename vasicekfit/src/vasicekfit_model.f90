! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Dmitriy Mayorov
module vasicekfit_model
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use vasicekfit_kinds, only : dp
   use vasicekfit_normal, only : normal_cdf, normal_quantile
   use vasicekfit_linalg, only : ols_regression, sample_variance
   implicit none
   private

   type, public :: vasicek_fit_result
      logical :: ok = .false.
      character(len=:), allocatable :: message
      integer :: n_observations = 0
      integer :: n_predictors = 0
      integer :: portfolio_size = 0
      logical :: bias_correct = .false.
      real(dp) :: p = 0.0_dp
      real(dp) :: rho = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: sigma2_raw = 0.0_dp
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: kappa(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: response_probit(:)
      real(dp), allocatable :: adjusted_response(:)
      real(dp), allocatable :: design(:,:)
      real(dp), allocatable :: beta_covariance(:,:)
   end type vasicek_fit_result

   type, public :: prediction_result
      logical :: ok = .false.
      character(len=:), allocatable :: message
      real(dp), allocatable :: values(:,:)
   end type prediction_result

   public :: fit_vasicek, predict_link, predict_response
   public :: predict_quantiles, coefficients

contains

   function fit_vasicek(response, predictors, bias_correct, portfolio_size) result(fit)
      real(dp), intent(in) :: response(:)
      real(dp), intent(in), optional :: predictors(:,:)
      logical, intent(in), optional :: bias_correct
      integer, intent(in), optional :: portfolio_size
      type(vasicek_fit_result) :: fit
      real(dp), allocatable :: y(:), x(:,:), beta(:), fitted(:), residuals(:), beta_covariance(:,:)
      real(dp) :: p0, vr, v0, denominator
      integer :: n, m, i
      logical :: regression_ok, use_bias

      fit%message = ''
      n = size(response)
      if (present(predictors)) then
         m = size(predictors,2)
         if (size(predictors,1) /= n) then
            fit%message = 'predictors and response have incompatible dimensions'
            return
         end if
      else
         m = 0
      end if
      if (n <= m + 1) then
         fit%message = 'more observations than regression coefficients are required'
         return
      end if
      if (.not. all(ieee_is_finite(response))) then
         fit%message = 'response values must be finite'
         return
      end if
      if (present(predictors)) then
         if (.not. all(ieee_is_finite(predictors))) then
            fit%message = 'predictor values must be finite'
            return
         end if
      end if

      y = response
      fit%portfolio_size = 0
      if (present(portfolio_size)) then
         if (portfolio_size < 2) then
            fit%message = 'portfolio_size must be at least 2'
            return
         end if
         if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            fit%message = 'response values must lie in [0, 1] with portfolio correction'
            return
         end if
         p0 = sum(y) / real(n,dp)
         vr = sample_variance(y)
         if (vr <= 0.0_dp) then
            fit%message = 'finite-portfolio correction requires positive sample variance'
            return
         end if
         v0 = vr - (p0 * (1.0_dp - p0) - vr) / real(portfolio_size - 1,dp)
         if (v0 <= 0.0_dp) then
            fit%message = 'finite-portfolio correction yielded non-positive variance'
            return
         end if
         y = p0 + (y - p0) * sqrt(v0 / vr)
         fit%portfolio_size = portfolio_size
      end if

      if (any(y <= 0.0_dp) .or. any(y >= 1.0_dp)) then
         fit%message = 'response values must lie strictly inside (0, 1)'
         return
      end if

      allocate(x(n,m+1))
      x(:,1) = 1.0_dp
      if (m > 0) x(:,2:m+1) = predictors
      fit%response_probit = normal_quantile(y)
      call ols_regression(x, fit%response_probit, beta, fitted, residuals, beta_covariance, regression_ok)
      if (.not. regression_ok) then
         fit%message = 'the probit-space regression design is singular'
         return
      end if

      fit%sigma2_raw = dot_product(residuals, residuals) / real(n,dp)
      fit%sigma2 = fit%sigma2_raw
      use_bias = .false.
      if (present(bias_correct)) use_bias = bias_correct
      if (use_bias) then
         if (n - m - 1 <= 0) then
            fit%message = 'bias correction has non-positive residual degrees of freedom'
            return
         end if
         fit%sigma2 = fit%sigma2 * real(n,dp) / real(n - m - 1,dp)
      end if

      denominator = sqrt(1.0_dp + fit%sigma2)
      fit%p = normal_cdf(beta(1) / denominator)
      fit%rho = fit%sigma2 / (1.0_dp + fit%sigma2)
      allocate(fit%kappa(m))
      if (m > 0) fit%kappa = beta(2:m+1) / denominator

      fit%n_observations = n
      fit%n_predictors = m
      fit%bias_correct = use_bias
      fit%beta = beta
      fit%fitted = fitted
      fit%residuals = residuals
      fit%adjusted_response = y
      fit%design = x
      fit%beta_covariance = beta_covariance
      fit%ok = .true.

      do i = 1, size(fit%kappa)
         if (.not. ieee_is_finite(fit%kappa(i))) then
            fit%ok = .false.
            fit%message = 'non-finite transformed parameter estimate'
            return
         end if
      end do
   end function fit_vasicek

   function coefficients(fit) result(values)
      type(vasicek_fit_result), intent(in) :: fit
      real(dp), allocatable :: values(:)
      allocate(values(fit%n_predictors + 2))
      values(1) = fit%p
      values(2) = fit%rho
      if (fit%n_predictors > 0) values(3:) = fit%kappa
   end function coefficients

   function predict_link(fit, predictors) result(prediction)
      type(vasicek_fit_result), intent(in) :: fit
      real(dp), intent(in), optional :: predictors(:,:)
      type(prediction_result) :: prediction
      real(dp), allocatable :: x(:,:)
      integer :: n

      prediction%message = ''
      if (.not. fit%ok) then
         prediction%message = 'cannot predict from an invalid fit'
         return
      end if
      if (present(predictors)) then
         n = size(predictors,1)
         if (size(predictors,2) /= fit%n_predictors) then
            prediction%message = 'new predictor matrix has the wrong number of columns'
            return
         end if
         allocate(x(n,fit%n_predictors+1))
         x(:,1) = 1.0_dp
         if (fit%n_predictors > 0) x(:,2:) = predictors
         allocate(prediction%values(n,1))
         prediction%values(:,1) = matmul(x, fit%beta)
      else
         allocate(prediction%values(fit%n_observations,1))
         prediction%values(:,1) = fit%fitted
      end if
      prediction%ok = .true.
   end function predict_link

   function predict_response(fit, predictors) result(prediction)
      type(vasicek_fit_result), intent(in) :: fit
      real(dp), intent(in), optional :: predictors(:,:)
      type(prediction_result) :: prediction
      real(dp), allocatable :: contribution(:)
      integer :: n

      prediction%message = ''
      if (.not. fit%ok) then
         prediction%message = 'cannot predict from an invalid fit'
         return
      end if
      if (present(predictors)) then
         n = size(predictors,1)
         if (size(predictors,2) /= fit%n_predictors) then
            prediction%message = 'new predictor matrix has the wrong number of columns'
            return
         end if
         allocate(contribution(n))
         if (fit%n_predictors > 0) then
            contribution = matmul(predictors, fit%kappa)
         else
            contribution = 0.0_dp
         end if
      else
         n = fit%n_observations
         allocate(contribution(n))
         if (fit%n_predictors > 0) then
            contribution = matmul(fit%design(:,2:), fit%kappa)
         else
            contribution = 0.0_dp
         end if
      end if
      allocate(prediction%values(n,1))
      prediction%values(:,1) = normal_cdf(normal_quantile(fit%p) + contribution)
      prediction%ok = .true.
   end function predict_response

   function predict_quantiles(fit, alpha, predictors, response_scale) result(prediction)
      type(vasicek_fit_result), intent(in) :: fit
      real(dp), intent(in) :: alpha(:)
      real(dp), intent(in), optional :: predictors(:,:)
      logical, intent(in), optional :: response_scale
      type(prediction_result) :: prediction
      real(dp), allocatable :: contribution(:)
      real(dp) :: base, sq_rho, sq_one_minus_rho
      integer :: n, i, j
      logical :: return_response

      prediction%message = ''
      if (.not. fit%ok) then
         prediction%message = 'cannot predict from an invalid fit'
         return
      end if
      if (size(alpha) == 0 .or. any(alpha <= 0.0_dp) .or. any(alpha >= 1.0_dp)) then
         prediction%message = 'alpha values must lie strictly inside (0, 1)'
         return
      end if
      if (present(predictors)) then
         n = size(predictors,1)
         if (size(predictors,2) /= fit%n_predictors) then
            prediction%message = 'new predictor matrix has the wrong number of columns'
            return
         end if
         allocate(contribution(n))
         if (fit%n_predictors > 0) then
            contribution = matmul(predictors, fit%kappa)
         else
            contribution = 0.0_dp
         end if
      else
         n = fit%n_observations
         allocate(contribution(n))
         if (fit%n_predictors > 0) then
            contribution = matmul(fit%design(:,2:), fit%kappa)
         else
            contribution = 0.0_dp
         end if
      end if

      return_response = .true.
      if (present(response_scale)) return_response = response_scale
      base = normal_quantile(fit%p)
      sq_rho = sqrt(fit%rho)
      sq_one_minus_rho = sqrt(1.0_dp - fit%rho)
      allocate(prediction%values(n,size(alpha)))
      do j = 1, size(alpha)
         do i = 1, n
            prediction%values(i,j) = (base + contribution(i) + sq_rho * normal_quantile(alpha(j))) / &
               sq_one_minus_rho
         end do
      end do
      if (return_response) prediction%values = normal_cdf(prediction%values)
      prediction%ok = .true.
   end function predict_quantiles

end module vasicekfit_model
