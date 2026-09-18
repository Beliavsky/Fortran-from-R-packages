! SPDX-License-Identifier: MIT
! SPDX-FileComment: Generalized linear models adapted from R's IRLS implementation.
module r_stats_glm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_distributions, only: r_pnorm, r_qnorm
   use r_kinds, only: dp
   use r_linalg, only: inverse_matrix, least_squares
   use r_stats_design, only: independent_columns
   use r_stats_types, only: glm_fit_t
   implicit none
   private

   integer, parameter, public :: glm_family_binomial = 1, glm_family_poisson = 2
   integer, parameter, public :: glm_family_gaussian = 3, glm_family_gamma = 4
   integer, parameter, public :: glm_family_quasibinomial = 5, glm_family_quasipoisson = 6
   integer, parameter, public :: glm_link_logit = 1, glm_link_log = 2
   integer, parameter, public :: glm_link_identity = 3, glm_link_inverse = 4
   integer, parameter, public :: glm_link_probit = 5, glm_link_cloglog = 6
   integer, parameter, public :: glm_link_cauchit = 7, glm_link_sqrt = 8
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt_two_pi = sqrt(2.0_dp*pi)

   public :: glm_binomial_fit, glm_fit, glm_gamma_fit, glm_gaussian_fit
   public :: glm_inverse_link, glm_link_value, glm_mu_eta
   public :: glm_pearson_resid, glm_poisson_fit, glm_predict_response
   public :: glm_quasibinomial_fit, glm_quasipoisson_fit

   interface glm_binomial_fit
      module procedure glm_binomial_fit_integer, glm_binomial_fit_real
   end interface glm_binomial_fit

   interface glm_poisson_fit
      module procedure glm_poisson_fit_integer, glm_poisson_fit_real
   end interface glm_poisson_fit

contains

   function glm_binomial_fit_integer(y, xpred, weights, offset, intercept, link) result(fit)
      !! Fits a binomial-logit GLM to integer zero-one responses.
      integer, intent(in) :: y(:) !! Binary response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Binomial link identifier; defaults to logit.
      type(glm_fit_t) :: fit !! Fitted binomial-logit model.
      integer :: selected_link

      selected_link = glm_link_logit
      if (present(link)) selected_link = link
      fit = fit_with_options(real(y, dp), xpred, glm_family_binomial, selected_link, &
                             weights, offset, intercept)
   end function glm_binomial_fit_integer

   function glm_binomial_fit_real(y, xpred, weights, offset, intercept, link) result(fit)
      !! Fits a binomial-logit GLM to proportions by iteratively reweighted least squares.
      real(dp), intent(in) :: y(:) !! Proportion response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: weights(:) !! Prior or binomial-trial weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Binomial link identifier; defaults to logit.
      type(glm_fit_t) :: fit !! Fitted binomial-logit model.
      integer :: selected_link

      selected_link = glm_link_logit
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, glm_family_binomial, selected_link, weights, offset, intercept)
   end function glm_binomial_fit_real

   function glm_poisson_fit_integer(y, xpred, offset, weights, intercept, link) result(fit)
      !! Fits a Poisson-log GLM to integer count responses.
      integer, intent(in) :: y(:) !! Nonnegative count response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Poisson link identifier; defaults to log.
      type(glm_fit_t) :: fit !! Fitted Poisson-log model.
      integer :: selected_link

      selected_link = glm_link_log
      if (present(link)) selected_link = link
      fit = fit_with_options(real(y, dp), xpred, glm_family_poisson, selected_link, &
                             weights, offset, intercept)
   end function glm_poisson_fit_integer

   function glm_poisson_fit_real(y, xpred, offset, weights, intercept, link) result(fit)
      !! Fits a Poisson-log GLM by iteratively reweighted least squares.
      real(dp), intent(in) :: y(:) !! Nonnegative response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Poisson link identifier; defaults to log.
      type(glm_fit_t) :: fit !! Fitted Poisson-log model.
      integer :: selected_link

      selected_link = glm_link_log
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, glm_family_poisson, selected_link, weights, offset, intercept)
   end function glm_poisson_fit_real

   function glm_quasibinomial_fit(y, xpred, weights, offset, intercept, link) result(fit)
      !! Fits a quasi-binomial GLM with estimated dispersion.
      real(dp), intent(in) :: y(:) !! Proportion response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Binomial link identifier; defaults to logit.
      type(glm_fit_t) :: fit !! Fitted quasi-binomial model.
      integer :: selected_link

      selected_link = glm_link_logit
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, glm_family_quasibinomial, selected_link, &
                             weights, offset, intercept)
   end function glm_quasibinomial_fit

   function glm_quasipoisson_fit(y, xpred, weights, offset, intercept, link) result(fit)
      !! Fits a quasi-Poisson GLM with estimated dispersion.
      real(dp), intent(in) :: y(:) !! Nonnegative response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Poisson link identifier; defaults to log.
      type(glm_fit_t) :: fit !! Fitted quasi-Poisson model.
      integer :: selected_link

      selected_link = glm_link_log
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, glm_family_quasipoisson, selected_link, &
                             weights, offset, intercept)
   end function glm_quasipoisson_fit

   function glm_gaussian_fit(y, xpred, weights, offset, intercept, link) result(fit)
      !! Fits a Gaussian GLM with identity, log, or inverse link.
      real(dp), intent(in) :: y(:) !! Real response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Link identifier; defaults to identity.
      type(glm_fit_t) :: fit !! Fitted Gaussian model.
      integer :: selected_link

      selected_link = glm_link_identity
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, glm_family_gaussian, selected_link, weights, offset, intercept)
   end function glm_gaussian_fit

   function glm_gamma_fit(y, xpred, weights, offset, intercept, link) result(fit)
      !! Fits a Gamma GLM with inverse, log, or identity link.
      real(dp), intent(in) :: y(:) !! Strictly positive response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      integer, intent(in), optional :: link !! Link identifier; defaults to inverse.
      type(glm_fit_t) :: fit !! Fitted Gamma model.
      integer :: selected_link

      selected_link = glm_link_inverse
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, glm_family_gamma, selected_link, weights, offset, intercept)
   end function glm_gamma_fit

   function glm_fit(y, xpred, family, link, weights, offset, intercept) result(fit)
      !! Fits a supported GLM family and link through the common typed API.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      integer, intent(in) :: family !! One of the public `glm_family_*` constants.
      integer, intent(in), optional :: link !! Link constant, or the family's default when absent.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      type(glm_fit_t) :: fit !! Fitted generalized linear model.
      integer :: selected_link

      selected_link = default_link(family)
      if (present(link)) selected_link = link
      fit = fit_with_options(y, xpred, family, selected_link, weights, offset, intercept)
   end function glm_fit

   function fit_with_options(y, xpred, family, link, weights, offset, intercept) result(fit)
      !! Materializes optional arguments before entering the numerical GLM engine.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      integer, intent(in) :: family !! Supported family identifier.
      integer, intent(in) :: link !! Supported link identifier.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in), optional :: offset(:) !! Additive linear-predictor offset with size `n`.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      type(glm_fit_t) :: fit !! Fitted generalized linear model.
      real(dp), allocatable :: actual_offset(:), prior_weights(:)
      logical :: has_intercept
      integer :: n

      n = size(y)
      if (size(xpred, 1) /= n) error stop "glm_fit: response and predictor row counts differ"
      if (.not. all(ieee_is_finite(xpred))) error stop "glm_fit: predictors must be finite"
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) error stop "glm_fit: invalid weights"
         prior_weights = weights
      else
         allocate (prior_weights(n), source=1.0_dp)
      end if
      if (present(offset)) then
         if (size(offset) /= n) error stop "glm_fit: offset size mismatch"
         if (.not. all(ieee_is_finite(offset))) error stop "glm_fit: offset must be finite"
         actual_offset = offset
      else
         allocate (actual_offset(n), source=0.0_dp)
      end if
      has_intercept = .true.
      if (present(intercept)) has_intercept = intercept
      fit = fit_glm_core(y, xpred, family, link, prior_weights, actual_offset, has_intercept)
      fit%weighted = present(weights)
   end function fit_with_options

   function fit_glm_core(y, xpred, family, link, prior_weights, offset, has_intercept) result(fit)
      !! Performs rank-aware iteratively reweighted least squares.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)`, excluding an intercept.
      integer, intent(in) :: family !! Supported family identifier.
      integer, intent(in) :: link !! Supported link identifier.
      real(dp), intent(in) :: prior_weights(:) !! Nonnegative prior weights with size `n`.
      real(dp), intent(in) :: offset(:) !! Additive offset with size `n`.
      logical, intent(in) :: has_intercept !! Whether to prepend an intercept column.
      type(glm_fit_t) :: fit !! Fitted generalized linear model.
      real(dp), allocatable :: beta(:), beta_new(:), design(:, :), design_rank(:, :)
      real(dp), allocatable :: eta(:), initial_eta(:), mu(:), response_work(:)
      real(dp), allocatable :: weighted_design(:, :), weighted_response(:), working_weights(:)
      real(dp) :: deviance_current, deviance_new
      integer, allocatable :: independent(:)
      integer :: half_step, info, iteration, j, k, n

      n = size(y)
      call validate_family_data(y, family, link, prior_weights)
      design = make_design(xpred, has_intercept)
      k = size(design, 2)
      if (k == 0) error stop "glm_fit: model has no coefficients"
      weighted_design = design
      do j = 1, n
         weighted_design(j, :) = sqrt(prior_weights(j))*weighted_design(j, :)
      end do
      call independent_columns(weighted_design, sqrt(epsilon(1.0_dp)), independent)
      fit%rank = size(independent)
      if (fit%rank == 0) error stop "glm_fit: zero-rank design matrix"
      if (count(prior_weights > 0.0_dp) <= fit%rank) error stop "glm_fit: no residual degrees of freedom"
      design_rank = design(:, independent)

      initial_eta = initial_linear_predictor(y, family, link) - offset
      weighted_design = design_rank
      weighted_response = initial_eta
      do j = 1, n
         weighted_design(j, :) = sqrt(prior_weights(j))*weighted_design(j, :)
         weighted_response(j) = sqrt(prior_weights(j))*weighted_response(j)
      end do
      allocate (beta(fit%rank), beta_new(fit%rank))
      call least_squares(weighted_design, weighted_response, beta, info)
      if (info /= 0) error stop "glm_fit: initial least-squares solution failed"

      allocate (eta(n), mu(n), response_work(n), working_weights(n))
      do iteration = 1, 100
         eta = offset + matmul(design_rank, beta)
         if (.not. valid_linear_predictor(eta, link)) exit
         mu = glm_inverse_link(eta, link)
         if (.not. valid_means(mu, family)) exit
         deviance_current = residual_deviance(y, mu, prior_weights, family)
         working_weights = prior_weights*glm_mu_eta(eta, link)**2/family_variance(mu, family)
         where (prior_weights > 0.0_dp)
            working_weights = max(working_weights, tiny(1.0_dp))
         elsewhere
            working_weights = 0.0_dp
         end where
         response_work = eta - offset + (y - mu)/glm_mu_eta(eta, link)
         weighted_design = design_rank
         weighted_response = response_work
         do j = 1, n
            weighted_design(j, :) = sqrt(working_weights(j))*weighted_design(j, :)
            weighted_response(j) = sqrt(working_weights(j))*weighted_response(j)
         end do
         call least_squares(weighted_design, weighted_response, beta_new, info)
         if (info /= 0) exit
         do half_step = 1, 25
            eta = offset + matmul(design_rank, beta_new)
            if (valid_linear_predictor(eta, link)) then
               mu = glm_inverse_link(eta, link)
               if (valid_means(mu, family)) exit
            end if
            beta_new = 0.5_dp*(beta + beta_new)
         end do
         if (.not. valid_linear_predictor(eta, link)) exit
         if (.not. valid_means(mu, family)) exit
         deviance_new = residual_deviance(y, mu, prior_weights, family)
         if (abs(deviance_new - deviance_current)/(0.1_dp + abs(deviance_new)) < 1.0e-8_dp .or. &
             maxval(abs(beta_new - beta)) < 1.0e-10_dp*(1.0_dp + maxval(abs(beta)))) then
            beta = beta_new
            fit%convergence = 0
            exit
         end if
         beta = beta_new
      end do

      fit%iter = min(iteration, 100)
      fit%family = family
      fit%link = link
      fit%has_intercept = has_intercept
      fit%weights = prior_weights
      fit%offset = offset
      fit%y = y
      fit%xpred = xpred
      allocate (fit%coef(k), source=0.0_dp)
      fit%coef(independent) = beta
      allocate (fit%aliased(k), source=.true.)
      fit%aliased(independent) = .false.
      call finish_glm(fit, design, independent)
   end function fit_glm_core

   pure function glm_predict_response(fit, xpred, offset) result(response)
      !! Predicts response means from a supported generalized linear model.
      type(glm_fit_t), intent(in) :: fit !! Previously fitted generalized linear model.
      real(dp), intent(in) :: xpred(:, :) !! New predictor rows, excluding an intercept column.
      real(dp), intent(in), optional :: offset(:) !! New additive offsets, one per predictor row.
      real(dp), allocatable :: response(:) !! Predicted response means.
      real(dp), allocatable :: eta(:)
      integer :: n, p

      n = size(xpred, 1)
      p = size(xpred, 2)
      if (fit%has_intercept) then
         if (size(fit%coef) /= p + 1) error stop "glm_predict_response: predictor count mismatch"
         eta = fit%coef(1) + matmul(xpred, fit%coef(2:))
      else
         if (size(fit%coef) /= p) error stop "glm_predict_response: predictor count mismatch"
         eta = matmul(xpred, fit%coef)
      end if
      if (present(offset)) then
         if (size(offset) /= n) error stop "glm_predict_response: offset size mismatch"
         eta = eta + offset
      else if (allocated(fit%offset)) then
         if (size(fit%offset) == n) eta = eta + fit%offset
      end if
      response = glm_inverse_link(eta, fit%link)
   end function glm_predict_response

   pure function glm_pearson_resid(fit) result(residuals)
      !! Computes prior-weighted Pearson residuals for a supported generalized linear model.
      type(glm_fit_t), intent(in) :: fit !! Previously fitted generalized linear model.
      real(dp), allocatable :: residuals(:) !! Pearson residuals with size `n`.

      residuals = sqrt(fit%weights)*(fit%y - fit%fitted)/ &
                  sqrt(family_variance(fit%fitted, fit%family))
   end function glm_pearson_resid

   subroutine finish_glm(fit, design, independent)
      !! Populates fitted values, residuals, covariance statistics, and diagnostics.
      type(glm_fit_t), intent(inout) :: fit !! Model populated in place after coefficient iteration.
      real(dp), intent(in) :: design(:, :) !! Full design matrix with shape `(n,p)`.
      integer, intent(in) :: independent(:) !! Indices of estimable design columns.
      real(dp), allocatable :: covariance_rank(:, :), design_rank(:, :), information(:, :)
      real(dp), allocatable :: eta(:), variance(:)
      real(dp) :: missing
      integer :: info, j, k

      k = size(design, 2)
      eta = fit%offset + matmul(design, fit%coef)
      fit%fitted = glm_inverse_link(eta, fit%link)
      fit%resid = fit%y - fit%fitted
      fit%df = count(fit%weights > 0.0_dp) - fit%rank
      variance = family_variance(fit%fitted, fit%family)
      fit%working_weights = fit%weights*glm_mu_eta(eta, fit%link)**2/variance
      fit%deviance = residual_deviance(fit%y, fit%fitted, fit%weights, fit%family)
      if (fit%family == glm_family_gaussian .or. fit%family == glm_family_gamma .or. &
          fit%family == glm_family_quasibinomial .or. fit%family == glm_family_quasipoisson) then
         fit%dispersion = sum(fit%weights*fit%resid**2/variance)/real(fit%df, dp)
      else
         fit%dispersion = 1.0_dp
      end if

      design_rank = design(:, independent)
      information = matmul(transpose(design_rank), &
                           spread(fit%working_weights, 2, fit%rank)*design_rank)
      call inverse_matrix(information, covariance_rank, info)
      missing = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate (fit%covariance(k, k), source=missing)
      allocate (fit%se(k), source=missing)
      allocate (fit%z_value(k), source=missing)
      allocate (fit%p_value(k), source=missing)
      if (info /= 0) return
      covariance_rank = fit%dispersion*covariance_rank
      do j = 1, fit%rank
         fit%covariance(independent, independent(j)) = covariance_rank(:, j)
         fit%se(independent(j)) = sqrt(max(0.0_dp, covariance_rank(j, j)))
      end do
      fit%z_value(independent) = fit%coef(independent)/max(fit%se(independent), tiny(1.0_dp))
      fit%p_value(independent) = 2.0_dp*r_pnorm(abs(fit%z_value(independent)), lower_tail=.false.)
   end subroutine finish_glm

   pure function make_design(xpred, has_intercept) result(design)
      !! Builds a design matrix with an optional leading intercept column.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix with shape `(n,p)`.
      logical, intent(in) :: has_intercept !! Whether to prepend a unit column.
      real(dp), allocatable :: design(:, :) !! Resulting design matrix.

      allocate (design(size(xpred, 1), size(xpred, 2) + merge(1, 0, has_intercept)))
      if (has_intercept) then
         design(:, 1) = 1.0_dp
         if (size(xpred, 2) > 0) design(:, 2:) = xpred
      else if (size(xpred, 2) > 0) then
         design = xpred
      end if
   end function make_design

   pure subroutine validate_family_data(y, family, link, weights)
      !! Rejects response, family, link, and weight combinations outside the supported domain.
      real(dp), intent(in) :: y(:) !! Response values to validate.
      integer, intent(in) :: family !! Family identifier to validate.
      integer, intent(in) :: link !! Link identifier to validate.
      real(dp), intent(in) :: weights(:) !! Prior weights to validate.

      if (size(weights) /= size(y) .or. any(weights < 0.0_dp)) error stop "glm_fit: invalid weights"
      if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(weights))) then
         error stop "glm_fit: response and weights must be finite"
      end if
      if (count(weights > 0.0_dp) == 0) error stop "glm_fit: no positive-weight observations"
      select case (family)
      case (glm_family_binomial, glm_family_quasibinomial)
         if (link /= glm_link_logit .and. link /= glm_link_probit .and. &
             link /= glm_link_cloglog .and. link /= glm_link_cauchit) then
            error stop "glm_fit: unsupported binomial link"
         end if
         if (any(y < 0.0_dp .or. y > 1.0_dp)) error stop "glm_fit: binomial responses outside [0,1]"
      case (glm_family_poisson, glm_family_quasipoisson)
         if (link /= glm_link_log .and. link /= glm_link_identity .and. link /= glm_link_sqrt) then
            error stop "glm_fit: unsupported Poisson link"
         end if
         if (any(y < 0.0_dp)) error stop "glm_fit: Poisson responses must be nonnegative"
      case (glm_family_gaussian)
         if (link /= glm_link_identity .and. link /= glm_link_log .and. link /= glm_link_inverse) then
            error stop "glm_fit: unsupported Gaussian link"
         end if
      case (glm_family_gamma)
         if (link /= glm_link_inverse .and. link /= glm_link_log .and. link /= glm_link_identity) then
            error stop "glm_fit: unsupported Gamma link"
         end if
         if (any(y <= 0.0_dp)) error stop "glm_fit: Gamma responses must be positive"
      case default
         error stop "glm_fit: unsupported family"
      end select
   end subroutine validate_family_data

   pure function initial_linear_predictor(y, family, link) result(eta)
      !! Constructs finite response-based starting values for IRLS.
      real(dp), intent(in) :: y(:) !! Validated response vector.
      integer, intent(in) :: family !! Supported family identifier.
      integer, intent(in) :: link !! Supported link identifier.
      real(dp), allocatable :: eta(:) !! Initial linear predictor with size `n`.
      real(dp), allocatable :: mu(:)

      select case (family)
      case (glm_family_binomial, glm_family_quasibinomial)
         mu = (y + 0.5_dp)/2.0_dp
      case (glm_family_poisson, glm_family_quasipoisson)
         mu = y + 0.1_dp
      case (glm_family_gamma)
         mu = y
      case default
         if (link == glm_link_log) then
            mu = max(y, sqrt(epsilon(1.0_dp)))
         else
            mu = y
         end if
      end select
      eta = glm_link_value(mu, link)
   end function initial_linear_predictor

   pure elemental function glm_inverse_link(eta, link) result(mu)
      !! Applies a supported inverse-link function with finite exponential bounds.
      real(dp), intent(in) :: eta !! Linear-predictor value.
      integer, intent(in) :: link !! Supported link identifier.
      real(dp) :: mu !! Mean-scale value.

      select case (link)
      case (glm_link_logit)
         mu = 1.0_dp/(1.0_dp + exp(-max(-35.0_dp, min(35.0_dp, eta))))
      case (glm_link_log)
         mu = exp(max(-35.0_dp, min(35.0_dp, eta)))
      case (glm_link_identity)
         mu = eta
      case (glm_link_inverse)
         mu = 1.0_dp/eta
      case (glm_link_probit)
         mu = r_pnorm(eta)
      case (glm_link_cloglog)
         mu = 1.0_dp - exp(-exp(max(-35.0_dp, min(35.0_dp, eta))))
      case (glm_link_cauchit)
         mu = 0.5_dp + atan(eta)/pi
      case (glm_link_sqrt)
         mu = eta**2
      case default
         mu = eta
      end select
   end function glm_inverse_link

   pure elemental function glm_link_value(mu, link) result(eta)
      !! Applies a supported link function.
      real(dp), intent(in) :: mu !! Mean-scale value in the link domain.
      integer, intent(in) :: link !! Supported link identifier.
      real(dp) :: eta !! Linear-predictor value.

      select case (link)
      case (glm_link_logit)
         eta = log(mu/(1.0_dp - mu))
      case (glm_link_log)
         eta = log(mu)
      case (glm_link_identity)
         eta = mu
      case (glm_link_inverse)
         eta = 1.0_dp/mu
      case (glm_link_probit)
         eta = r_qnorm(mu)
      case (glm_link_cloglog)
         eta = log(-log(1.0_dp - mu))
      case (glm_link_cauchit)
         eta = tan(pi*(mu - 0.5_dp))
      case (glm_link_sqrt)
         eta = sqrt(mu)
      case default
         eta = mu
      end select
   end function glm_link_value

   pure elemental function glm_mu_eta(eta, link) result(derivative)
      !! Returns the derivative of the mean with respect to the linear predictor.
      real(dp), intent(in) :: eta !! Linear-predictor value.
      integer, intent(in) :: link !! Supported link identifier.
      real(dp) :: derivative !! Value of `d mu / d eta`.
      real(dp) :: mu

      mu = glm_inverse_link(eta, link)

      select case (link)
      case (glm_link_logit)
         derivative = max(tiny(1.0_dp), mu*(1.0_dp - mu))
      case (glm_link_log)
         derivative = max(tiny(1.0_dp), mu)
      case (glm_link_identity)
         derivative = 1.0_dp
      case (glm_link_inverse)
         derivative = -mu**2
      case (glm_link_probit)
         derivative = max(tiny(1.0_dp), exp(-0.5_dp*eta**2)/sqrt_two_pi)
      case (glm_link_cloglog)
         derivative = max(tiny(1.0_dp), exp(max(-35.0_dp, min(35.0_dp, eta)) - &
                                            exp(max(-35.0_dp, min(35.0_dp, eta)))))
      case (glm_link_cauchit)
         derivative = max(tiny(1.0_dp), 1.0_dp/(pi*(1.0_dp + eta**2)))
      case (glm_link_sqrt)
         derivative = 2.0_dp*eta
      case default
         derivative = 1.0_dp
      end select
   end function glm_mu_eta

   pure elemental function family_variance(mu, family) result(variance)
      !! Evaluates the supported GLM variance function.
      real(dp), intent(in) :: mu !! Mean-scale value.
      integer, intent(in) :: family !! Supported family identifier.
      real(dp) :: variance !! Unit-dispersion variance value.

      select case (family)
      case (glm_family_binomial, glm_family_quasibinomial)
         variance = max(tiny(1.0_dp), mu*(1.0_dp - mu))
      case (glm_family_poisson, glm_family_quasipoisson)
         variance = max(tiny(1.0_dp), mu)
      case (glm_family_gamma)
         variance = max(tiny(1.0_dp), mu**2)
      case default
         variance = 1.0_dp
      end select
   end function family_variance

   pure function residual_deviance(y, mu, weights, family) result(value)
      !! Computes the family-specific residual deviance.
      real(dp), intent(in) :: y(:) !! Observed responses with size `n`.
      real(dp), intent(in) :: mu(:) !! Fitted means with size `n`.
      real(dp), intent(in) :: weights(:) !! Prior weights with size `n`.
      integer, intent(in) :: family !! Supported family identifier.
      real(dp) :: value !! Sum of weighted unit-deviance contributions.
      real(dp) :: contribution
      integer :: i

      value = 0.0_dp
      do i = 1, size(y)
         select case (family)
         case (glm_family_binomial, glm_family_quasibinomial)
            contribution = 0.0_dp
            if (y(i) > 0.0_dp) contribution = contribution + y(i)*log(y(i)/mu(i))
            if (y(i) < 1.0_dp) then
               contribution = contribution + (1.0_dp - y(i))*log((1.0_dp - y(i))/(1.0_dp - mu(i)))
            end if
            value = value + 2.0_dp*weights(i)*contribution
         case (glm_family_poisson, glm_family_quasipoisson)
            if (y(i) > 0.0_dp) then
               contribution = y(i)*log(y(i)/mu(i)) - (y(i) - mu(i))
            else
               contribution = mu(i)
            end if
            value = value + 2.0_dp*weights(i)*contribution
         case (glm_family_gamma)
            value = value + 2.0_dp*weights(i)*((y(i) - mu(i))/mu(i) - log(y(i)/mu(i)))
         case default
            value = value + weights(i)*(y(i) - mu(i))**2
         end select
      end do
   end function residual_deviance

   pure function valid_means(mu, family) result(valid)
      !! Reports whether every fitted mean lies in the selected family's domain.
      real(dp), intent(in) :: mu(:) !! Candidate fitted means.
      integer, intent(in) :: family !! Supported family identifier.
      logical :: valid !! True when all means satisfy the family domain.

      select case (family)
      case (glm_family_binomial, glm_family_quasibinomial)
         valid = all(ieee_is_finite(mu) .and. mu > 0.0_dp .and. mu < 1.0_dp)
      case (glm_family_poisson, glm_family_gamma, glm_family_quasipoisson)
         valid = all(ieee_is_finite(mu) .and. mu > 0.0_dp)
      case default
         valid = all(ieee_is_finite(mu))
      end select
   end function valid_means

   pure function valid_linear_predictor(eta, link) result(valid)
      !! Reports whether a linear predictor lies in the selected inverse-link domain.
      real(dp), intent(in) :: eta(:) !! Candidate linear-predictor values.
      integer, intent(in) :: link !! Supported link identifier.
      logical :: valid !! True when the inverse link is finite and defined.

      valid = all(ieee_is_finite(eta))
      if (link == glm_link_inverse) valid = valid .and. all(abs(eta) > tiny(1.0_dp))
      if (link == glm_link_sqrt) valid = valid .and. all(abs(eta) > tiny(1.0_dp))
   end function valid_linear_predictor

   pure function default_link(family) result(link)
      !! Returns R's default link identifier for a supported family.
      integer, intent(in) :: family !! Supported family identifier.
      integer :: link !! Default link identifier.

      select case (family)
      case (glm_family_binomial, glm_family_quasibinomial)
         link = glm_link_logit
      case (glm_family_poisson, glm_family_quasipoisson)
         link = glm_link_log
      case (glm_family_gaussian)
         link = glm_link_identity
      case (glm_family_gamma)
         link = glm_link_inverse
      case default
         error stop "glm_fit: unsupported family"
      end select
   end function default_link

end module r_stats_glm
