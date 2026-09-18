! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_core
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_control, flexmix_result
   use flexmix_types, only : flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial
   use flexmix_types, only : flexmix_model_mvnorm, flexmix_model_mvbinary, flexmix_model_mvpois, flexmix_model_mvcombi
   use flexmix_types, only : flexmix_model_lognormal, flexmix_model_exponential, flexmix_model_inverse_gaussian
   use flexmix_types, only : flexmix_model_gamma, flexmix_model_weibull
   use flexmix_types, only : flexmix_model_gamma_regression, flexmix_model_multinomial
   use flexmix_types, only : flexmix_model_ziglm_poisson, flexmix_model_ziglm_binomial
   use flexmix_types, only : flexmix_model_robust_gaussian, flexmix_model_robust_poisson, flexmix_model_factanal
   use flexmix_types, only : flexmix_model_conditional_logit
   use flexmix_types, only : flexmix_role_regular, flexmix_role_structural_zero, flexmix_role_robust_background
   use flexmix_components, only : fit_gaussian_regression_component, gaussian_regression_log_density
   use flexmix_components, only : fit_poisson_regression_component, poisson_regression_log_density
   use flexmix_components, only : fit_binomial_regression_component, binomial_regression_log_density
   use flexmix_components, only : fit_gamma_regression_component, gamma_regression_log_density
   use flexmix_components, only : fit_multinomial_regression_component, multinomial_regression_log_density
   use flexmix_components, only : fit_conditional_logit_component, conditional_logit_log_density
   use flexmix_components, only : fit_mvnormal_component, fit_factor_analysis_component, mvnormal_log_density
   use flexmix_components, only : fit_mvbinary_component, mvbinary_log_density
   use flexmix_components, only : fit_mvpois_component, mvpois_log_density
   use flexmix_components, only : fit_mvcombi_component, mvcombi_log_density
   use flexmix_components, only : fit_lognormal_component, lognormal_log_density
   use flexmix_components, only : fit_exponential_component, exponential_log_density
   use flexmix_components, only : fit_inverse_gaussian_component, inverse_gaussian_log_density
   use flexmix_components, only : fit_gamma_component, gamma_log_density
   use flexmix_components, only : fit_weibull_component, weibull_log_density
   use flexmix_penalized, only : fit_glmnet_component
   use flexmix_smooth, only : fit_smooth_component
   use flexmix_em_utils, only : initialize_posteriors, classify_posteriors, group_first_mask, group_log_densities
   use flexmix_em_utils, only : fit_constant_prior, fit_multinomial_prior, multinomial_prior
   use flexmix_em_utils, only : select_active_components, weighted_loglikelihood, argmax_rows
   use flexmix_numeric, only : normalize_log_probabilities, weighted_variance_unbiased
   implicit none
   private
   public :: flexmix_fit_regression
   public :: flexmix_fit_glmnet
   public :: flexmix_fit_mgcv
   public :: flexmix_fit_binomial
   public :: flexmix_fit_multinomial
   public :: flexmix_fit_conditional_logit
   public :: flexmix_fit_ziglm_poisson
   public :: flexmix_fit_ziglm_binomial
   public :: flexmix_fit_robust_gaussian
   public :: flexmix_fit_robust_poisson
   public :: flexmix_fit_multivariate
   public :: flexmix_fit_univariate

contains

   subroutine flexmix_fit_regression(x, y, k, model_kind, result, control, case_weights, group, concomitant_x, &
                                     initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`, including any desired intercept column.
      real(dp), intent(in) :: y(:) !! Response vector, size `n`; Gaussian values or nonnegative Poisson counts.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      integer, intent(in) :: model_kind !! Regression family identifier: Gaussian or Poisson.
      type(flexmix_result), intent(out) :: result !! Fitted mixture result with active components packed contiguously.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior matrix `(n, q)`, constant within groups.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), sigma(:), gamma_shape(:), work(:), old_beta(:)
      real(dp) :: llh, old_llh
      logical :: changed
      integer :: component_df, fit_info, info, iter, j, n, p

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x, 2)
      call initialize_result_header(result, model_kind, n, p, 1, k)
      if (size(x, 1) /= n .or. k < 1) then
         result%status = -1
         return
      end if
      if (model_kind /= flexmix_model_gaussian .and. model_kind /= flexmix_model_poisson .and. &
          model_kind /= flexmix_model_gamma_regression) then
         result%status = -2
         return
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), beta(p,k), sigma(k), &
               gamma_shape(k), work(n), old_beta(p))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = -5
            return
         end if
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      beta = 0.0_dp
      sigma = 1.0_dp
      gamma_shape = 1.0_dp
      prior = 1.0_dp / real(k, dp)
      prior_rows = spread(prior, 1, n)
      if (present(concomitant_x)) then
         if (size(concomitant_x, 1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2), k))
         coef_prior = 0.0_dp
      else
         allocate(coef_prior(1,k))
         coef_prior = 0.0_dp
      end if
      call initialize_posteriors(y, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1, ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            old_beta = beta(:,j)
            select case (model_kind)
            case (flexmix_model_gaussian)
               call fit_gaussian_regression_component(x, y, work, beta(:,j), sigma(j), component_df, fit_info, offset)
            case (flexmix_model_poisson)
               call fit_poisson_regression_component(x, y, work, beta(:,j), component_df, fit_info, old_beta, offset)
            case (flexmix_model_gamma_regression)
               call fit_gamma_regression_component(x, y, work, beta(:,j), gamma_shape(j), component_df, fit_info, old_beta, offset)
            end select
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         call evaluate_regression_components(x, y, model_kind, active, beta, sigma, gamma_shape, log_density, offset)
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter, max(1, ctl%iter_max))
      result%loglik = llh
      call pack_regression_result(result, active, posterior, logp, beta, sigma, weights, first, &
                                  coef_prior, present(concomitant_x), gamma_shape)
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_regression

   subroutine flexmix_fit_glmnet(x, y, k, model_kind, result, trials, control, case_weights, group, concomitant_x, &
                                 initial_cluster, initial_posterior, offset, adaptive, select, alpha, nfolds, fold_id, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)` whose first column is an intercept of ones.
      real(dp), intent(in) :: y(:) !! Gaussian/Poisson response, or binomial success counts for a binomial model.
      integer, intent(in) :: k !! Requested initial number of penalized mixture components, at least one.
      integer, intent(in) :: model_kind !! Family identifier: Gaussian, Poisson, or binomial.
      type(flexmix_result), intent(out) :: result !! Fitted penalized-regression mixture with selected penalties by retained.
      real(dp), intent(in), optional :: trials(:) !! Binomial trial totals, size `n`; required only for binomial models.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix `(n,q)`, constant within.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      logical, intent(in), optional :: adaptive !! Use adaptive-lasso penalty factors from unpenalized component fits; default.
      logical, intent(in), optional :: select(:) !! Penalized-term mask for columns two through `p`, size `p-1`; false terms.
      real(dp), intent(in), optional :: alpha !! Elastic-net mixing fraction in `[0,1]`; default one (lasso).
      integer, intent(in), optional :: nfolds !! Number of fixed cross-validation folds when `fold_id` is absent; default five.
      integer, intent(in), optional :: fold_id(:) !! Optional fixed one-based cross-validation fold labels, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative cross-validation penalty grid.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), sigma(:), work(:), signal(:), work_trials(:), selected_lambda(:), packed_lambda(:)
      integer, allocatable :: component_df(:)
      real(dp) :: alpha_value, llh, old_llh
      logical :: changed
      integer :: a, fit_info, info, iter, j, n, p, concomitant_df

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      call initialize_result_header(result, model_kind, n, p, 1, k)
      if (size(x,1) /= n .or. k < 1) then
         result%status = -1
         return
      end if
      if (model_kind /= flexmix_model_gaussian .and. model_kind /= flexmix_model_poisson .and. &
          model_kind /= flexmix_model_binomial) then
         result%status = -2
         return
      end if
      if (model_kind == flexmix_model_binomial) then
         if (.not. present(trials)) then
            result%status = -5
            return
         end if
         if (size(trials) /= n .or. any(trials < y) .or. any(y < 0.0_dp)) then
            result%status = -6
            return
         end if
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = -7
            return
         end if
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), beta(p,k), sigma(k), work(n))
      allocate(signal(n), work_trials(n), selected_lambda(k), component_df(k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      beta = 0.0_dp
      sigma = 1.0_dp
      selected_lambda = 0.0_dp
      component_df = 0
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2),k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      work_trials = 1.0_dp
      signal = y
      if (model_kind == flexmix_model_binomial) then
         work_trials = trials
         signal = y / max(work_trials,1.0_dp)
      end if
      call initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified,active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            call fit_glmnet_component(x, y, trials, work, model_kind, beta(:,j), sigma(j), component_df(j), &
                                      selected_lambda(j), fit_info, adaptive, select, alpha, nfolds, fold_id, offset, lambda_grid)
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active,prior)
         do j = 1, k
            if (.not. active(j)) then
               log_density(:,j) = -huge(1.0_dp)
            else
               select case (model_kind)
               case (flexmix_model_gaussian)
                  call gaussian_regression_log_density(x, y, beta(:,j), sigma(j), log_density(:,j), offset)
               case (flexmix_model_poisson)
                  call poisson_regression_log_density(x, y, beta(:,j), log_density(:,j), offset)
               case (flexmix_model_binomial)
                  call binomial_regression_log_density(x, y, work_trials, beta(:,j), log_density(:,j), offset)
               end select
            end if
         end do
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active,prior,concomitant_x,coef_prior,prior_rows)
         call combine_log_density(grouped_density,prior_rows,active,logp)
         call normalize_log_probabilities(logp,posterior,logsum)
         llh = weighted_loglikelihood(logsum,weights,first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      call pack_regression_result(result, active, posterior, logp, beta, sigma, weights, first, coef_prior, present(concomitant_x))
      allocate(packed_lambda(result%k))
      a = 0
      result%df = 0
      do j = 1, k
         if (.not. active(j)) cycle
         a = a + 1
         packed_lambda(a) = selected_lambda(j)
         result%df = result%df + component_df(j)
      end do
      allocate(result%penalty_lambda(result%k))
      result%penalty_lambda = packed_lambda
      alpha_value = 1.0_dp
      if (present(alpha)) alpha_value = min(1.0_dp,max(0.0_dp,alpha))
      result%penalty_alpha = alpha_value
      if (present(concomitant_x)) then
         concomitant_df = size(coef_prior,1) * max(0,result%k - 1)
      else
         concomitant_df = max(0,result%k - 1)
      end if
      result%df = result%df + concomitant_df
      if (model_kind /= flexmix_model_gaussian .and. allocated(result%sigma)) deallocate(result%sigma)
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_glmnet

   subroutine flexmix_fit_mgcv(x, y, k, model_kind, penalty_matrix, result, trials, control, case_weights, group, &
                               concomitant_x, initial_cluster, initial_posterior, offset, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Prefit GAM/design matrix `(n,p)` containing parametric and smooth-basis columns.
      real(dp), intent(in) :: y(:) !! Gaussian/Poisson response, or binomial success counts for a binomial family.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      integer, intent(in) :: model_kind !! Family identifier: Gaussian, Poisson, or binomial.
      real(dp), intent(in) :: penalty_matrix(:,:) !! Symmetric smoothing penalty matrix `(p,p)` aligned with `x` columns.
      type(flexmix_result), intent(out) :: result !! Fitted penalized-smooth mixture with per-component smoothing diagnostics.
      real(dp), intent(in), optional :: trials(:) !! Binomial trial totals, size `n`; required only for binomial models.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix `(n,q)`, constant within.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative smoothing-parameter grid.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), sigma(:), work(:), signal(:), work_trials(:), selected_lambda(:), edf(:)
      real(dp), allocatable :: packed_lambda(:), packed_edf(:)
      real(dp) :: llh, old_llh
      logical :: changed
      integer :: a, fit_info, info, iter, j, n, p, concomitant_df

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      call initialize_result_header(result, model_kind, n, p, 1, k)
      if (size(x,1) /= n .or. k < 1 .or. any(shape(penalty_matrix) /= [p,p])) then
         result%status = -1
         return
      end if
      if (model_kind /= flexmix_model_gaussian .and. model_kind /= flexmix_model_poisson .and. &
          model_kind /= flexmix_model_binomial) then
         result%status = -2
         return
      end if
      if (model_kind == flexmix_model_binomial) then
         if (.not. present(trials)) then
            result%status = -5
            return
         end if
         if (size(trials) /= n .or. any(trials < y) .or. any(y < 0.0_dp)) then
            result%status = -6
            return
         end if
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = -7
            return
         end if
      end if
      allocate(weights(n),first(n),active(k),posterior(n,k),classified(n,k),prior(k),prior_rows(n,k))
      allocate(log_density(n,k),grouped_density(n,k),logp(n,k),logsum(n),beta(p,k),sigma(k),work(n))
      allocate(signal(n),work_trials(n),selected_lambda(k),edf(k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      call make_first_mask(n,group,first,info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      beta = 0.0_dp
      sigma = 1.0_dp
      selected_lambda = 0.0_dp
      edf = 0.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2),k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      work_trials = 1.0_dp
      signal = y
      if (model_kind == flexmix_model_binomial) then
         work_trials = trials
         signal = y / max(work_trials,1.0_dp)
      end if
      call initialize_posteriors(signal,k,posterior,initial_cluster,initial_posterior,info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior,ctl%classify,classified)
         call update_prior(classified,weights,first,active,ctl%minprior,concomitant_x,coef_prior, &
                           prior,prior_rows,changed,info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified,active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            call fit_smooth_component(x, y, trials, work, model_kind, penalty_matrix, beta(:,j), sigma(j), edf(j), &
                                      selected_lambda(j), fit_info, lambda_grid, offset)
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active,prior)
         do j = 1, k
            if (.not. active(j)) then
               log_density(:,j) = -huge(1.0_dp)
            else
               select case (model_kind)
               case (flexmix_model_gaussian)
                  call gaussian_regression_log_density(x,y,beta(:,j),sigma(j),log_density(:,j),offset)
               case (flexmix_model_poisson)
                  call poisson_regression_log_density(x,y,beta(:,j),log_density(:,j),offset)
               case (flexmix_model_binomial)
                  call binomial_regression_log_density(x,y,work_trials,beta(:,j),log_density(:,j),offset)
               end select
            end if
         end do
         call apply_grouping(log_density,group,grouped_density)
         call refresh_prior_rows(active,prior,concomitant_x,coef_prior,prior_rows)
         call combine_log_density(grouped_density,prior_rows,active,logp)
         call normalize_log_probabilities(logp,posterior,logsum)
         llh = weighted_loglikelihood(logsum,weights,first)
         if (iter > 1) then
            if (abs(llh-old_llh)/(abs(llh)+0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      call pack_regression_result(result,active,posterior,logp,beta,sigma,weights,first,coef_prior,present(concomitant_x))
      allocate(packed_lambda(result%k),packed_edf(result%k))
      a = 0
      result%effective_df = 0.0_dp
      do j = 1, k
         if (.not. active(j)) cycle
         a = a + 1
         packed_lambda(a) = selected_lambda(j)
         packed_edf(a) = edf(j)
         result%effective_df = result%effective_df + edf(j)
      end do
      allocate(result%smoothing_lambda(result%k),result%component_effective_df(result%k))
      result%smoothing_lambda = packed_lambda
      result%component_effective_df = packed_edf
      if (present(concomitant_x)) then
         concomitant_df = size(coef_prior,1) * max(0,result%k-1)
      else
         concomitant_df = max(0,result%k-1)
      end if
      result%effective_df = result%effective_df + real(concomitant_df,dp)
      result%df = nint(result%effective_df)
      if (model_kind /= flexmix_model_gaussian .and. allocated(result%sigma)) deallocate(result%sigma)
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_mgcv

   subroutine flexmix_fit_binomial(x, success, failure, k, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`, including any desired intercept column.
      real(dp), intent(in) :: success(:) !! Nonnegative binomial success counts, size `n`.
      real(dp), intent(in) :: failure(:) !! Nonnegative binomial failure counts, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      type(flexmix_result), intent(out) :: result !! Fitted binomial-regression mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive logit-scale offset, size `n`.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), work(:), old_beta(:), trials(:), signal(:)
      real(dp) :: llh, old_llh
      logical :: changed
      integer :: component_df, fit_info, info, iter, j, n, p

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(success)
      p = size(x, 2)
      call initialize_result_header(result, flexmix_model_binomial, n, p, 1, k)
      if (size(failure) /= n .or. size(x,1) /= n .or. k < 1 .or. any(success < 0.0_dp) .or. any(failure < 0.0_dp)) then
         result%status = -1
         return
      end if
      allocate(trials(n), signal(n))
      trials = success + failure
      signal = success / max(trials, 1.0_dp)
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), beta(p,k), work(n), old_beta(p))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = -5
            return
         end if
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      beta = 0.0_dp
      prior = 1.0_dp / real(k, dp)
      prior_rows = spread(prior, 1, n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2), k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      call initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1, ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            old_beta = beta(:,j)
            call fit_binomial_regression_component(x, success, trials, work, beta(:,j), component_df, fit_info, old_beta, offset)
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (active(j)) then
               call binomial_regression_log_density(x, success, trials, beta(:,j), log_density(:,j), offset)
            else
               log_density(:,j) = -huge(1.0_dp)
            end if
         end do
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter, max(1, ctl%iter_max))
      result%loglik = llh
      call pack_regression_result(result, active, posterior, logp, beta, [(1.0_dp, j = 1, k)], &
                                  weights, first, coef_prior, present(concomitant_x))
      if (allocated(result%sigma)) deallocate(result%sigma)
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_binomial

   subroutine flexmix_fit_multinomial(x, y_class, k, result, control, case_weights, group, concomitant_x, &
                                      initial_cluster, initial_posterior)
      real(dp), intent(in) :: x(:,:) !! Multinomial-regression design matrix `(n, p)`, with an intercept column if wanted.
      integer, intent(in) :: y_class(:) !! One-based class labels, size `n`, covering every integer through `max(y_class)`.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      type(flexmix_result), intent(out) :: result !! Fitted mixture of multinomial-logit regressions.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting mixture-component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting mixture posterior probabilities, shape `(n.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:), class_seen(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: coef(:,:,:), old_coef(:,:), work(:), signal(:)
      real(dp) :: llh, old_llh
      logical :: changed
      integer :: component_df, fit_info, info, iter, j, n, p, nclass, c

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y_class)
      p = size(x,2)
      if (n > 0) then
         nclass = maxval(y_class)
      else
         nclass = 0
      end if
      call initialize_result_header(result, flexmix_model_multinomial, n, p, nclass, k)
      result%nclass = nclass
      if (size(x,1) /= n .or. k < 1 .or. nclass < 2 .or. any(y_class < 1)) then
         result%status = -1
         return
      end if
      allocate(class_seen(nclass))
      class_seen = .false.
      do j = 1, n
         if (y_class(j) <= nclass) class_seen(y_class(j)) = .true.
      end do
      if (.not. all(class_seen)) then
         result%status = -2
         return
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), work(n), signal(n))
      allocate(coef(p,nclass-1,k), old_coef(p,nclass-1))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      coef = 0.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2),k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      signal = real(y_class,dp)
      call initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            old_coef = coef(:,:,j)
            call fit_multinomial_regression_component(x, y_class, nclass, work, coef(:,:,j), component_df, fit_info, old_coef)
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (active(j)) then
               call multinomial_regression_log_density(x, y_class, coef(:,:,j), log_density(:,j))
            else
               log_density(:,j) = -huge(1.0_dp)
            end if
         end do
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      call pack_multinomial_result(result, active, posterior, logp, coef, weights, first, coef_prior, present(concomitant_x))
      result%nclass = nclass
      c = count(active)
      if (c < 1) result%status = 2
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_multinomial

   subroutine flexmix_fit_conditional_logit(x, y, strata, k, result, control, case_weights, &
                                               initial_cluster, initial_posterior)
      real(dp), intent(in) :: x(:,:) !! Conditional-logit covariate matrix without an intercept, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Binary choice indicator, size `n`; exactly one selected alternative is required per.
      integer, intent(in) :: strata(:) !! Integer choice-set labels, size `n`; equal values identify one conditional likelihood.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      type(flexmix_result), intent(out) :: result !! Fitted conditional-logit mixture with coefficients and posteriors.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, classification, and minimum-prior.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative stratum-frequency weights, constant within each.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based labels, size `n`; preferably constant by group.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n,k)`.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), old_beta(:), work(:), signal(:), sigma_dummy(:), group_mean(:)
      real(dp) :: llh, old_llh, tolerance
      logical :: changed
      integer :: component_df, fit_info, info, iter, j, n, p, i, r, group_size, event_count

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      call initialize_result_header(result, flexmix_model_conditional_logit, n, p, 1, k)
      if (size(x,1) /= n .or. size(strata) /= n .or. k < 1 .or. p < 1) then
         result%status = -1
         return
      end if
      if (any(abs(y) > 100.0_dp * epsilon(1.0_dp) .and. &
              abs(y - 1.0_dp) > 100.0_dp * epsilon(1.0_dp))) then
         result%status = -2
         return
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), work(n), signal(n))
      allocate(beta(p,k), old_beta(p), sigma_dummy(k), group_mean(k), coef_prior(1,k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      call group_first_mask(strata, first)
      do i = 1, n
         if (.not. first(i)) cycle
         group_size = 0
         event_count = 0
         do r = 1, n
            if (strata(r) /= strata(i)) cycle
            group_size = group_size + 1
            if (y(r) > 0.5_dp) event_count = event_count + 1
            tolerance = 1.0e-10_dp * (1.0_dp + max(abs(weights(i)),abs(weights(r))))
            if (abs(weights(r) - weights(i)) > tolerance) then
               result%status = -4
               return
            end if
         end do
         if (group_size < 2 .or. event_count /= 1) then
            result%status = -5
            return
         end if
      end do
      active = .true.
      beta = 0.0_dp
      sigma_dummy = 1.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      coef_prior = 0.0_dp
      signal = real(strata,dp)
      call initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      do i = 1, n
         if (.not. first(i)) cycle
         group_mean = 0.0_dp
         group_size = 0
         do r = 1, n
            if (strata(r) /= strata(i)) cycle
            group_mean = group_mean + posterior(r,:)
            group_size = group_size + 1
         end do
         group_mean = group_mean / real(group_size,dp)
         do r = 1, n
            if (strata(r) == strata(i)) posterior(r,:) = group_mean
         end do
      end do
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, coef=coef_prior, prior=prior, &
                           prior_rows=prior_rows, changed=changed, info=info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            old_beta = beta(:,j)
            call fit_conditional_logit_component(x, y, strata, work, beta(:,j), component_df, fit_info, old_beta)
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (active(j)) then
               call conditional_logit_log_density(x, y, strata, beta(:,j), log_density(:,j), fit_info)
               if (fit_info /= 0) active(j) = .false.
            else
               log_density(:,j) = -huge(1.0_dp)
            end if
         end do
         call ensure_one_active(active, prior)
         call group_log_densities(log_density, strata, grouped_density)
         prior_rows = spread(prior,1,n)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      call pack_regression_result(result, active, posterior, logp, beta, sigma_dummy, weights, first, coef_prior, .false.)
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_conditional_logit

   subroutine flexmix_fit_ziglm_poisson(x, y, k, result, control, case_weights, group, concomitant_x, &
                                       initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Poisson regression design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson count response, size `n`.
      integer, intent(in) :: k !! Requested initial number of components including the structural-zero component.
      type(flexmix_result), intent(out) :: result !! Fitted zero-inflated Poisson mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive log-link offset, size `n`.
      call flexmix_fit_special_scalar(x, y, k, flexmix_model_ziglm_poisson, result, control, case_weights, group, &
                                      concomitant_x, initial_cluster, initial_posterior, offset, .false.)
   end subroutine flexmix_fit_ziglm_poisson

   subroutine flexmix_fit_robust_gaussian(x, y, k, result, bgw, control, case_weights, group, concomitant_x, &
                                          initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Gaussian regression design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      integer, intent(in) :: k !! Requested initial number of components including the robust background component.
      type(flexmix_result), intent(out) :: result !! Fitted robust Gaussian mixture-of-regressions result.
      logical, intent(in), optional :: bgw !! Use current posterior weights for the background fit when true; default false.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive identity-link offset, size `n`.
      logical :: use_bgw
      use_bgw = .false.
      if (present(bgw)) use_bgw = bgw
      call flexmix_fit_special_scalar(x, y, k, flexmix_model_robust_gaussian, result, control, case_weights, group, &
                                      concomitant_x, initial_cluster, initial_posterior, offset, use_bgw)
   end subroutine flexmix_fit_robust_gaussian

   subroutine flexmix_fit_robust_poisson(x, y, k, result, bgw, control, case_weights, group, concomitant_x, &
                                         initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Poisson regression design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson count response, size `n`.
      integer, intent(in) :: k !! Requested initial number of components including the robust background component.
      type(flexmix_result), intent(out) :: result !! Fitted robust Poisson mixture-of-regressions result.
      logical, intent(in), optional :: bgw !! Use current posterior weights for the background fit when true; default false.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive log-link offset, size `n`.
      logical :: use_bgw
      use_bgw = .false.
      if (present(bgw)) use_bgw = bgw
      call flexmix_fit_special_scalar(x, y, k, flexmix_model_robust_poisson, result, control, case_weights, group, &
                                      concomitant_x, initial_cluster, initial_posterior, offset, use_bgw)
   end subroutine flexmix_fit_robust_poisson

   subroutine flexmix_fit_special_scalar(x, y, k, special_kind, result, control, case_weights, group, concomitant_x, &
                                         initial_cluster, initial_posterior, offset, bgw)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Scalar response vector, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      integer, intent(in) :: special_kind !! Specialized family identifier for zero-inflated or robust Gaussian/Poisson fitting.
      type(flexmix_result), intent(out) :: result !! Fitted specialized regression mixture.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive family-link offset, size `n`.
      logical, intent(in) :: bgw !! Whether the robust background fit uses its current posterior weights.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), sigma(:), work(:), old_beta(:), background_weights(:)
      real(dp) :: llh, old_llh, background_mean, background_var, sw
      logical :: changed, variance_ok, special_retained
      integer :: base_kind, component_df, fit_info, info, iter, j, n, p, concomitant_df, ordinary_df, a

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      select case (special_kind)
      case (flexmix_model_ziglm_poisson, flexmix_model_robust_poisson)
         base_kind = flexmix_model_poisson
      case (flexmix_model_robust_gaussian)
         base_kind = flexmix_model_gaussian
      case default
         call initialize_result_header(result, special_kind, n, p, 1, k)
         result%status = -2
         return
      end select
      call initialize_result_header(result, special_kind, n, p, 1, k)
      if (size(x,1) /= n .or. k < 1 .or. p < 1) then
         result%status = -1
         return
      end if
      if (any(abs(x(:,1) - 1.0_dp) > 100.0_dp * epsilon(1.0_dp))) then
         result%status = -6
         return
      end if
      if (base_kind == flexmix_model_poisson .and. any(y < 0.0_dp)) then
         result%status = -7
         return
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), beta(p,k), sigma(k), work(n), old_beta(p))
      allocate(background_weights(n))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = -5
            return
         end if
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      beta = 0.0_dp
      sigma = 1.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2),k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      call initialize_posteriors(y, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            if (j == 1) then
               select case (special_kind)
               case (flexmix_model_ziglm_poisson)
                  beta(:,j) = 0.0_dp
                  fit_info = 0
               case (flexmix_model_robust_gaussian)
                  if (bgw) then
                     background_weights = weights * classified(:,j)
                  else
                     background_weights = 1.0_dp
                  end if
                  call weighted_variance_unbiased(y, background_weights, background_mean, background_var, variance_ok)
                  if (variance_ok) then
                     beta(:,j) = 0.0_dp
                     beta(1,j) = background_mean
                     sigma(j) = sqrt(max(background_var,tiny(1.0_dp)))
                     fit_info = 0
                  else
                     fit_info = 1
                  end if
               case (flexmix_model_robust_poisson)
                  if (bgw) then
                     background_weights = weights * classified(:,j)
                  else
                     background_weights = 1.0_dp
                  end if
                  sw = sum(background_weights)
                  if (sw > tiny(1.0_dp)) then
                     background_mean = dot_product(background_weights,y) / sw
                     beta(:,j) = 0.0_dp
                     beta(1,j) = log(max(3.0_dp * background_mean,tiny(1.0_dp)))
                     fit_info = 0
                  else
                     fit_info = 1
                  end if
               end select
            else
               work = weights * classified(:,j)
               old_beta = beta(:,j)
               if (base_kind == flexmix_model_gaussian) then
                  call fit_gaussian_regression_component(x, y, work, beta(:,j), sigma(j), component_df, fit_info, offset)
               else
                  call fit_poisson_regression_component(x, y, work, beta(:,j), component_df, fit_info, old_beta, offset)
               end if
            end if
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (.not. active(j)) then
               log_density(:,j) = -huge(1.0_dp)
            else if (j == 1 .and. special_kind == flexmix_model_ziglm_poisson) then
               where (abs(y) <= 100.0_dp * epsilon(1.0_dp))
                  log_density(:,j) = 0.0_dp
               elsewhere
                  log_density(:,j) = -huge(1.0_dp)
               end where
            else if (base_kind == flexmix_model_gaussian) then
               call gaussian_regression_log_density(x, y, beta(:,j), sigma(j), log_density(:,j), offset)
            else
               call poisson_regression_log_density(x, y, beta(:,j), log_density(:,j), offset)
            end if
         end do
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      special_retained = active(1)
      result%model_kind = base_kind
      call pack_regression_result(result, active, posterior, logp, beta, sigma, weights, first, coef_prior, present(concomitant_x))
      allocate(result%component_role(result%k))
      result%component_role = flexmix_role_regular
      a = 0
      do j = 1, k
         if (.not. active(j)) cycle
         a = a + 1
         if (j == 1) then
            if (special_kind == flexmix_model_ziglm_poisson) then
               result%component_role(a) = flexmix_role_structural_zero
            else
               result%component_role(a) = flexmix_role_robust_background
            end if
         end if
      end do
      if (special_retained) result%model_kind = special_kind
      if (present(concomitant_x)) then
         concomitant_df = size(coef_prior,1) * max(0,result%k - 1)
      else
         concomitant_df = max(0,result%k - 1)
      end if
      if (base_kind == flexmix_model_gaussian) then
         ordinary_df = p + 1
      else
         ordinary_df = p
      end if
      if (special_retained) then
         result%df = max(0,result%k - 1) * ordinary_df + concomitant_df
      else
         result%df = result%k * ordinary_df + concomitant_df
      end if
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_special_scalar

   subroutine flexmix_fit_ziglm_binomial(x, success, failure, k, result, control, case_weights, group, concomitant_x, &
                                         initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Binomial regression design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: success(:) !! Nonnegative binomial success counts, size `n`.
      real(dp), intent(in) :: failure(:) !! Nonnegative binomial failure counts, size `n`.
      integer, intent(in) :: k !! Requested initial number of components including the structural-zero component.
      type(flexmix_result), intent(out) :: result !! Fitted zero-inflated binomial mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive logit-link offset, size `n`.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:)
      real(dp), allocatable :: coef_prior(:,:), log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: beta(:,:), work(:), old_beta(:), trials(:), signal(:), sigma_dummy(:)
      real(dp) :: llh, old_llh
      logical :: changed, special_retained
      integer :: component_df, fit_info, info, iter, j, n, p, concomitant_df, a

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(success)
      p = size(x,2)
      call initialize_result_header(result, flexmix_model_ziglm_binomial, n, p, 1, k)
      if (size(failure) /= n .or. size(x,1) /= n .or. k < 1 .or. p < 1 .or. &
          any(success < 0.0_dp) .or. any(failure < 0.0_dp)) then
         result%status = -1
         return
      end if
      if (any(abs(x(:,1) - 1.0_dp) > 100.0_dp * epsilon(1.0_dp))) then
         result%status = -6
         return
      end if
      allocate(trials(n), signal(n))
      trials = success + failure
      signal = success / max(trials,1.0_dp)
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), beta(p,k), work(n), old_beta(p), sigma_dummy(k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = -5
            return
         end if
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      beta = 0.0_dp
      sigma_dummy = 1.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2),k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      call initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            if (j == 1) then
               beta(:,j) = 0.0_dp
               fit_info = 0
            else
               work = weights * classified(:,j)
               old_beta = beta(:,j)
               call fit_binomial_regression_component(x, success, trials, work, beta(:,j), component_df, fit_info, old_beta, offset)
            end if
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (.not. active(j)) then
               log_density(:,j) = -huge(1.0_dp)
            else if (j == 1) then
               where (abs(success) <= 100.0_dp * epsilon(1.0_dp))
                  log_density(:,j) = 0.0_dp
               elsewhere
                  log_density(:,j) = -huge(1.0_dp)
               end where
            else
               call binomial_regression_log_density(x, success, trials, beta(:,j), log_density(:,j), offset)
            end if
         end do
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      special_retained = active(1)
      result%model_kind = flexmix_model_binomial
      call pack_regression_result(result, active, posterior, logp, beta, sigma_dummy, weights, first, coef_prior, &
                                  present(concomitant_x))
      allocate(result%component_role(result%k))
      result%component_role = flexmix_role_regular
      a = 0
      do j = 1, k
         if (.not. active(j)) cycle
         a = a + 1
         if (j == 1) result%component_role(a) = flexmix_role_structural_zero
      end do
      if (special_retained) result%model_kind = flexmix_model_ziglm_binomial
      if (present(concomitant_x)) then
         concomitant_df = size(coef_prior,1) * max(0,result%k - 1)
      else
         concomitant_df = max(0,result%k - 1)
      end if
      if (special_retained) then
         result%df = max(0,result%k - 1) * p + concomitant_df
      else
         result%df = result%k * p + concomitant_df
      end if
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_ziglm_binomial

   subroutine flexmix_fit_multivariate(y, k, model_kind, result, control, case_weights, group, concomitant_x, &
                                       initial_cluster, initial_posterior, diagonal, truncated, binary, factors)
      real(dp), intent(in) :: y(:,:) !! Multivariate response matrix, shape `(n, d)`.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      integer, intent(in) :: model_kind !! Component family: multivariate normal, binary, Poisson, or combined.
      type(flexmix_result), intent(out) :: result !! Fitted multivariate mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n, k)`.
      logical, intent(in), optional :: diagonal !! Use diagonal multivariate-normal covariance when true; default false.
      logical, intent(in), optional :: truncated !! For binary components, condition on at least one success when true.
      logical, intent(in), optional :: binary(:) !! Bernoulli-column mask; remaining columns are independent Gaussian.
      integer, intent(in), optional :: factors !! Number of latent factors for `flexmix_model_factanal`; ignored by other.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:), binary_mask(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:), coef_prior(:,:)
      real(dp), allocatable :: log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:), work(:), signal(:)
      real(dp), allocatable :: center(:,:), covariance(:,:,:), probability(:,:), lambda_value(:,:), variance(:,:)
      real(dp), allocatable :: factor_loadings(:,:,:), uniqueness(:,:), marginal_variance(:,:)
      real(dp) :: llh, old_llh
      logical :: changed, use_diagonal, use_truncated
      integer :: component_df, fit_info, info, iter, j, n, d, factor_count

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y,1)
      d = size(y,2)
      call initialize_result_header(result, model_kind, n, 0, d, k)
      if (n < 1 .or. d < 1 .or. k < 1) then
         result%status = -1
         return
      end if
      if ((model_kind < flexmix_model_mvnorm .or. model_kind > flexmix_model_mvcombi) .and. &
          model_kind /= flexmix_model_factanal) then
         result%status = -2
         return
      end if
      factor_count = 1
      if (model_kind == flexmix_model_factanal) then
         if (.not. present(factors)) then
            result%status = -7
            return
         end if
         factor_count = factors
         if (factor_count < 1 .or. factor_count >= d .or. (d-factor_count)**2-d-factor_count < 0) then
            result%status = -7
            return
         end if
      end if
      use_diagonal = .false.
      if (present(diagonal)) use_diagonal = diagonal
      use_truncated = .false.
      if (present(truncated)) use_truncated = truncated
      allocate(binary_mask(d))
      binary_mask = .false.
      if (present(binary)) then
         if (size(binary) /= d) then
            result%status = -5
            return
         end if
         binary_mask = binary
      end if
      if (model_kind == flexmix_model_mvcombi .and. .not. present(binary)) then
         result%status = -6
         return
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), work(n), signal(n))
      allocate(center(d,k), covariance(d,d,k), probability(d,k), lambda_value(d,k), variance(d,k))
      allocate(factor_loadings(d,factor_count,k), uniqueness(d,k), marginal_variance(d,k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      center = 0.0_dp
      covariance = 0.0_dp
      probability = 0.5_dp
      lambda_value = 1.0_dp
      variance = 1.0_dp
      factor_loadings = 0.0_dp
      uniqueness = 1.0_dp
      marginal_variance = 1.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      signal = y(:,1)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2), k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      call initialize_posteriors(signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1, ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            select case (model_kind)
            case (flexmix_model_mvnorm)
               call fit_mvnormal_component(y, work, center(:,j), covariance(:,:,j), use_diagonal, component_df, fit_info)
            case (flexmix_model_mvbinary)
               call fit_mvbinary_component(y, work, probability(:,j), use_truncated, component_df, fit_info)
            case (flexmix_model_mvpois)
               call fit_mvpois_component(y, work, lambda_value(:,j), component_df, fit_info)
            case (flexmix_model_mvcombi)
               call fit_mvcombi_component(y, work, binary_mask, center(:,j), variance(:,j), component_df, fit_info)
            case (flexmix_model_factanal)
               call fit_factor_analysis_component(y, work, factor_count, center(:,j), covariance(:,:,j), &
                                                  marginal_variance(:,j), factor_loadings(:,:,j), uniqueness(:,j), &
                                                  component_df, fit_info)
            end select
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (.not. active(j)) then
               log_density(:,j) = -huge(1.0_dp)
               cycle
            end if
            select case (model_kind)
            case (flexmix_model_mvnorm)
               call mvnormal_log_density(y, center(:,j), covariance(:,:,j), log_density(:,j), fit_info)
               if (fit_info /= 0) active(j) = .false.
            case (flexmix_model_mvbinary)
               call mvbinary_log_density(y, probability(:,j), use_truncated, log_density(:,j))
            case (flexmix_model_mvpois)
               call mvpois_log_density(y, lambda_value(:,j), log_density(:,j))
            case (flexmix_model_mvcombi)
               call mvcombi_log_density(y, binary_mask, center(:,j), variance(:,j), log_density(:,j))
            case (flexmix_model_factanal)
               call mvnormal_log_density(y, center(:,j), covariance(:,:,j), log_density(:,j), fit_info)
               if (fit_info /= 0) active(j) = .false.
            end select
         end do
         call ensure_one_active(active, prior)
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter, max(1,ctl%iter_max))
      result%loglik = llh
      call pack_multivariate_result(result, active, posterior, logp, center, covariance, probability, lambda_value, variance, &
                                    factor_loadings, uniqueness, marginal_variance, factor_count, weights, first, coef_prior, &
                                    present(concomitant_x), use_diagonal, use_truncated, binary_mask)
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_multivariate

   subroutine flexmix_fit_univariate(y, k, model_kind, result, control, case_weights, group, concomitant_x, &
                                     initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:) !! Positive response values for the selected univariate distribution, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components, at least one.
      integer, intent(in) :: model_kind !! Distribution identifier: lognormal, exponential, inverse Gaussian, Gamma, or Weibull.
      type(flexmix_result), intent(out) :: result !! Fitted univariate mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`; rows in one group share posterior.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax-prior design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior probabilities, shape `(n, k)`.
      type(flexmix_control) :: ctl
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), prior_rows(:,:), coef_prior(:,:)
      real(dp), allocatable :: log_density(:,:), grouped_density(:,:), logp(:,:), logsum(:), work(:)
      real(dp), allocatable :: par1(:), par2(:)
      real(dp) :: llh, old_llh
      logical :: changed
      integer :: component_df, fit_info, info, iter, j, n

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      call initialize_result_header(result, model_kind, n, 0, 1, k)
      if (n < 1 .or. k < 1 .or. any(y <= 0.0_dp)) then
         result%status = -1
         return
      end if
      if (model_kind < flexmix_model_lognormal .or. model_kind > flexmix_model_weibull) then
         result%status = -2
         return
      end if
      allocate(weights(n), first(n), active(k), posterior(n,k), classified(n,k), prior(k), prior_rows(n,k))
      allocate(log_density(n,k), grouped_density(n,k), logp(n,k), logsum(n), work(n), par1(k), par2(k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -3
            return
         end if
         weights = case_weights
      end if
      call make_first_mask(n, group, first, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      active = .true.
      par1 = 1.0_dp
      par2 = 1.0_dp
      prior = 1.0_dp / real(k,dp)
      prior_rows = spread(prior,1,n)
      if (present(concomitant_x)) then
         if (size(concomitant_x,1) /= n) then
            result%status = -4
            return
         end if
         allocate(coef_prior(size(concomitant_x,2), k))
      else
         allocate(coef_prior(1,k))
      end if
      coef_prior = 0.0_dp
      call initialize_posteriors(y, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1, ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         call update_prior(classified, weights, first, active, ctl%minprior, concomitant_x, coef_prior, &
                           prior, prior_rows, changed, info)
         if (info /= 0) then
            result%status = 20 + info
            exit
         end if
         call renormalize_active(classified, active)
         do j = 1, k
            if (.not. active(j)) cycle
            work = weights * classified(:,j)
            select case (model_kind)
            case (flexmix_model_lognormal)
               call fit_lognormal_component(y, work, par1(j), par2(j), component_df, fit_info)
            case (flexmix_model_exponential)
               call fit_exponential_component(y, work, par1(j), component_df, fit_info)
               par2(j) = 1.0_dp
            case (flexmix_model_inverse_gaussian)
               call fit_inverse_gaussian_component(y, work, par1(j), par2(j), component_df, fit_info)
            case (flexmix_model_gamma)
               call fit_gamma_component(y, work, par1(j), par2(j), component_df, fit_info)
            case (flexmix_model_weibull)
               call fit_weibull_component(y, work, par1(j), par2(j), component_df, fit_info)
            end select
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_one_active(active, prior)
         do j = 1, k
            if (.not. active(j)) then
               log_density(:,j) = -huge(1.0_dp)
               cycle
            end if
            select case (model_kind)
            case (flexmix_model_lognormal)
               call lognormal_log_density(y, par1(j), par2(j), log_density(:,j))
            case (flexmix_model_exponential)
               call exponential_log_density(y, par1(j), log_density(:,j))
            case (flexmix_model_inverse_gaussian)
               call inverse_gaussian_log_density(y, par1(j), par2(j), log_density(:,j))
            case (flexmix_model_gamma)
               call gamma_log_density(y, par1(j), par2(j), log_density(:,j))
            case (flexmix_model_weibull)
               call weibull_log_density(y, par1(j), par2(j), log_density(:,j))
            end select
         end do
         call apply_grouping(log_density, group, grouped_density)
         call refresh_prior_rows(active, prior, concomitant_x, coef_prior, prior_rows)
         call combine_log_density(grouped_density, prior_rows, active, logp)
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = weighted_loglikelihood(logsum, weights, first)
         if (iter > 1) then
            if (abs(llh - old_llh) / (abs(llh) + 0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter, max(1,ctl%iter_max))
      result%loglik = llh
      call pack_univariate_result(result, active, posterior, logp, par1, par2, weights, first, &
                                  coef_prior, present(concomitant_x))
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_univariate

   subroutine initialize_result_header(result, model_kind, n, p, d, k)
      type(flexmix_result), intent(out) :: result !! Result whose scalar metadata is initialized; arrays remain unallocated.
      integer, intent(in) :: model_kind !! Integer identifier of the requested component model family.
      integer, intent(in) :: n !! Number of observation rows in the fit.
      integer, intent(in) :: p !! Number of component regression coefficients, or zero for clustering models.
      integer, intent(in) :: d !! Response dimension for clustering/distribution models.
      integer, intent(in) :: k !! Initial requested number of mixture components.
      result%model_kind = model_kind
      result%n = n
      result%p = p
      result%d = d
      result%k0 = k
      result%k = 0
      result%df = 0
      result%iterations = 0
      result%status = 0
      result%converged = .false.
      result%loglik = -huge(1.0_dp)
   end subroutine initialize_result_header

   subroutine make_first_mask(n, group, first, info)
      integer, intent(in) :: n !! Number of rows for which the group-first mask is constructed.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      logical, intent(out) :: first(:) !! Output mask selecting the first row of each group, size `n`.
      integer, intent(out) :: info !! Zero on success; `-7` when a supplied group vector has the wrong length.
      if (present(group)) then
         if (size(group) /= n) then
            first = .false.
            info = -7
            return
         end if
         call group_first_mask(group, first)
      else
         first = .true.
      end if
      info = 0
   end subroutine make_first_mask

   subroutine update_prior(classified, weights, first, active, minprior, concomitant_x, coef, prior, prior_rows, changed, info)
      real(dp), intent(in) :: classified(:,:) !! M-step membership weights, shape `(n, k0)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask, size `n`.
      logical, intent(inout) :: active(:) !! Active-component mask, size `k0`, possibly reduced in place.
      real(dp), intent(in) :: minprior !! Minimum marginal component prior retained by the EM fit.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional softmax prior design matrix, shape `(n, q)`.
      real(dp), intent(inout) :: coef(:,:) !! Concomitant coefficients, shape `(q, k0)` when used.
      real(dp), intent(out) :: prior(:) !! Marginal component prior probabilities, size `k0`.
      real(dp), intent(out) :: prior_rows(:,:) !! Observation-specific prior probabilities, shape `(n, k0)`.
      logical, intent(out) :: changed !! True when one or more components were deleted by `minprior`.
      integer, intent(out) :: info !! Zero on success; positive value for concomitant Newton failure.
      integer :: i
      real(dp) :: sw
      info = 0
      if (present(concomitant_x)) then
         call fit_multinomial_prior(concomitant_x, classified, weights, first, active, coef, prior_rows, info)
         if (info /= 0) return
         prior = 0.0_dp
         sw = 0.0_dp
         do i = 1, size(first)
            if (first(i)) then
               prior = prior + weights(i) * prior_rows(i,:)
               sw = sw + weights(i)
            end if
         end do
         if (sw > 0.0_dp) prior = prior / sw
      else
         call fit_constant_prior(classified, weights, first, active, prior)
         prior_rows = spread(prior, 1, size(classified,1))
      end if
      call select_active_components(prior, minprior, active, changed)
      call normalize_prior(active, prior)
      if (present(concomitant_x)) then
         call multinomial_prior(concomitant_x, coef, active, prior_rows)
      else
         prior_rows = spread(prior, 1, size(classified,1))
      end if
   end subroutine update_prior

   pure subroutine normalize_prior(active, prior)
      logical, intent(in) :: active(:) !! Active-component mask identifying which prior entries participate in normalization.
      real(dp), intent(inout) :: prior(:) !! Marginal component priors normalized in place over active components.
      real(dp) :: total
      integer :: j
      do j = 1, size(prior)
         if (.not. active(j)) prior(j) = 0.0_dp
      end do
      total = sum(prior)
      if (total > 0.0_dp) then
         prior = prior / total
      else
         do j = 1, size(prior)
            if (active(j)) prior(j) = 1.0_dp / real(count(active), dp)
         end do
      end if
   end subroutine normalize_prior

   pure subroutine ensure_one_active(active, prior)
      logical, intent(inout) :: active(:) !! Component activity mask that is guaranteed to retain at least one entry.
      real(dp), intent(in) :: prior(:) !! Most recent component priors used to choose a fallback component.
      integer :: best
      if (count(active) == 0) then
         best = maxloc(prior, dim=1)
         active(best) = .true.
      end if
   end subroutine ensure_one_active

   pure subroutine renormalize_active(classified, active)
      real(dp), intent(inout) :: classified(:,:) !! Membership matrix renormalized over active components row by row.
      logical, intent(in) :: active(:) !! Active-component mask, size equal to the number of matrix columns.
      real(dp) :: total
      integer :: i, j
      do j = 1, size(active)
         if (.not. active(j)) classified(:,j) = 0.0_dp
      end do
      do i = 1, size(classified,1)
         total = sum(classified(i,:))
         if (total > 0.0_dp) then
            classified(i,:) = classified(i,:) / total
         else
            do j = 1, size(active)
               if (active(j)) classified(i,j) = 1.0_dp / real(count(active),dp)
            end do
         end if
      end do
   end subroutine renormalize_active

   subroutine refresh_prior_rows(active, prior, concomitant_x, coef, prior_rows)
      logical, intent(in) :: active(:) !! Active-component mask, size `k0`.
      real(dp), intent(in) :: prior(:) !! Constant marginal priors used when no concomitant model is supplied.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant design matrix, shape `(n, q)`.
      real(dp), intent(in) :: coef(:,:) !! Concomitant softmax coefficients, shape `(q, k0)` when used.
      real(dp), intent(out) :: prior_rows(:,:) !! Row-wise component priors, shape `(n, k0)`.
      if (present(concomitant_x)) then
         call multinomial_prior(concomitant_x, coef, active, prior_rows)
      else
         prior_rows = spread(prior, 1, size(prior_rows,1))
      end if
   end subroutine refresh_prior_rows

   subroutine evaluate_regression_components(x, y, model_kind, active, beta, sigma, shape_value, log_density, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Regression response vector, size `n`.
      integer, intent(in) :: model_kind !! Gaussian, Poisson, or Gamma regression family identifier.
      logical, intent(in) :: active(:) !! Active-component mask, size `k0`.
      real(dp), intent(in) :: beta(:,:) !! Component coefficient matrix, shape `(p, k0)`.
      real(dp), intent(in) :: sigma(:) !! Gaussian component scales, size `k0`; ignored by other families.
      real(dp), intent(in) :: shape_value(:) !! Gamma shape parameters, size `k0`; ignored by other families.
      real(dp), intent(out) :: log_density(:,:) !! Per-row component log densities, shape `(n, k0)`.
      real(dp), intent(in), optional :: offset(:) !! Optional family-scale linear-predictor offset, size `n`.
      integer :: j
      do j = 1, size(active)
         if (.not. active(j)) then
            log_density(:,j) = -huge(1.0_dp)
         else
            select case (model_kind)
            case (flexmix_model_gaussian)
               call gaussian_regression_log_density(x, y, beta(:,j), sigma(j), log_density(:,j), offset)
            case (flexmix_model_poisson)
               call poisson_regression_log_density(x, y, beta(:,j), log_density(:,j), offset)
            case (flexmix_model_gamma_regression)
               call gamma_regression_log_density(x, y, beta(:,j), shape_value(j), log_density(:,j), offset)
            end select
         end if
      end do
   end subroutine evaluate_regression_components

   subroutine apply_grouping(log_density, group, grouped_density)
      real(dp), intent(in) :: log_density(:,:) !! Observation-level component log densities, shape `(n, k0)`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(out) :: grouped_density(:,:) !! Group-summed densities repeated by group, shape `(n, k0)`.
      if (present(group)) then
         call group_log_densities(log_density, group, grouped_density)
      else
         grouped_density = log_density
      end if
   end subroutine apply_grouping

   pure subroutine combine_log_density(log_density, prior_rows, active, logp)
      real(dp), intent(in) :: log_density(:,:) !! Component log densities, shape `(n, k0)`.
      real(dp), intent(in) :: prior_rows(:,:) !! Row-wise component prior probabilities, shape `(n, k0)`.
      logical, intent(in) :: active(:) !! Active-component mask, size `k0`.
      real(dp), intent(out) :: logp(:,:) !! Unnormalized posterior log weights, shape `(n, k0)`.
      integer :: j
      do j = 1, size(active)
         if (active(j)) then
            logp(:,j) = log_density(:,j) + log(max(prior_rows(:,j), tiny(1.0_dp)))
         else
            logp(:,j) = -huge(1.0_dp)
         end if
      end do
   end subroutine combine_log_density

   subroutine pack_common(result, active, posterior, logp, weights, first, coef_prior, has_concomitant)
      type(flexmix_result), intent(inout) :: result !! Fit receiving packed posterior, prior, cluster, and concomitant fields.
      logical, intent(in) :: active(:) !! Active-component mask over the original `k0` components.
      real(dp), intent(in) :: posterior(:,:) !! Final posterior matrix before packing, shape `(n, k0)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized posterior log weights, shape `(n, k0)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask used for marginal priors and likelihood accounting.
      real(dp), intent(in) :: coef_prior(:,:) !! Original concomitant coefficient matrix, possibly a one-row placeholder.
      logical, intent(in) :: has_concomitant !! True when concomitant coefficients should be retained in the result.
      integer, allocatable :: map(:)
      real(dp) :: sw, max_log
      integer :: a, i, j, kk
      kk = count(active)
      result%k = kk
      allocate(map(kk))
      a = 0
      do j = 1, size(active)
         if (active(j)) then
            a = a + 1
            map(a) = j
         end if
      end do
      allocate(result%posterior(size(posterior,1),kk), result%posterior_unscaled(size(posterior,1),kk))
      allocate(result%log_posterior_unscaled(size(posterior,1),kk), result%prior(kk))
      allocate(result%cluster(size(posterior,1)), result%size(kk))
      allocate(result%group_first(size(first)), result%case_weights(size(weights)))
      result%group_first = first
      result%case_weights = weights
      do a = 1, kk
         result%posterior(:,a) = posterior(:,map(a))
         result%log_posterior_unscaled(:,a) = logp(:,map(a))
         max_log = log(huge(1.0_dp)) - 2.0_dp
         result%posterior_unscaled(:,a) = exp(max(-max_log, min(max_log, logp(:,map(a)))))
      end do
      do i = 1, size(result%posterior,1)
         if (sum(result%posterior(i,:)) > 0.0_dp) result%posterior(i,:) = result%posterior(i,:) / sum(result%posterior(i,:))
      end do
      result%prior = 0.0_dp
      sw = 0.0_dp
      do i = 1, size(first)
         if (first(i)) then
            result%prior = result%prior + weights(i) * result%posterior(i,:)
            sw = sw + weights(i)
         end if
      end do
      if (sw > 0.0_dp) result%prior = result%prior / sw
      call argmax_rows(result%posterior, result%cluster)
      result%size = 0
      do i = 1, size(result%cluster)
         result%size(result%cluster(i)) = result%size(result%cluster(i)) + nint(weights(i))
      end do
      if (has_concomitant) then
         allocate(result%concomitant_coef(size(coef_prior,1),kk))
         do a = 1, kk
            result%concomitant_coef(:,a) = coef_prior(:,map(a))
         end do
      end if
   end subroutine pack_common

   subroutine pack_regression_result(result, active, posterior, logp, beta, sigma, weights, first, coef_prior, has_concomitant, &
                                     shape_value)
      type(flexmix_result), intent(inout) :: result !! Regression fit result receiving packed parameters and common mixture.
      logical, intent(in) :: active(:) !! Active-component mask over original components.
      real(dp), intent(in) :: posterior(:,:) !! Final posterior probabilities, shape `(n, k0)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized posterior log weights, shape `(n, k0)`.
      real(dp), intent(in) :: beta(:,:) !! Component regression coefficients, shape `(p, k0)`.
      real(dp), intent(in) :: sigma(:) !! Component Gaussian scales or placeholder values, size `k0`.
      real(dp), intent(in) :: weights(:) !! Case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask, size `n`.
      real(dp), intent(in) :: coef_prior(:,:) !! Concomitant coefficient matrix or placeholder.
      logical, intent(in) :: has_concomitant !! True when concomitant coefficients should be stored.
      real(dp), intent(in), optional :: shape_value(:) !! Gamma regression shape values by original component, size `k0`.
      integer :: a, j, kk, component_df, concomitant_df
      call pack_common(result, active, posterior, logp, weights, first, coef_prior, has_concomitant)
      kk = result%k
      allocate(result%beta(size(beta,1),kk))
      if (result%model_kind == flexmix_model_gaussian) allocate(result%sigma(kk))
      if (result%model_kind == flexmix_model_gamma_regression) allocate(result%shape(kk))
      a = 0
      do j = 1, size(active)
         if (active(j)) then
            a = a + 1
            result%beta(:,a) = beta(:,j)
            if (allocated(result%sigma)) result%sigma(a) = sigma(j)
            if (allocated(result%shape) .and. present(shape_value)) result%shape(a) = shape_value(j)
         end if
      end do
      select case (result%model_kind)
      case (flexmix_model_gaussian, flexmix_model_gamma_regression)
         component_df = size(beta,1) + 1
      case default
         component_df = size(beta,1)
      end select
      if (has_concomitant) then
         concomitant_df = size(coef_prior,1) * max(0,kk - 1)
      else
         concomitant_df = max(0,kk - 1)
      end if
      result%df = kk * component_df + concomitant_df
   end subroutine pack_regression_result

   subroutine pack_multinomial_result(result, active, posterior, logp, coef, weights, first, coef_prior, has_concomitant)
      type(flexmix_result), intent(inout) :: result !! Multinomial mixture result receiving packed component parameters.
      logical, intent(in) :: active(:) !! Active-component mask over original mixture components.
      real(dp), intent(in) :: posterior(:,:) !! Final posterior probabilities, shape `(n, k0)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized posterior log weights, shape `(n, k0)`.
      real(dp), intent(in) :: coef(:,:,:) !! Multinomial coefficient tensor, shape `(p, nclass-1, k0)`.
      real(dp), intent(in) :: weights(:) !! Case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask, size `n`.
      real(dp), intent(in) :: coef_prior(:,:) !! Concomitant coefficient matrix or placeholder.
      logical, intent(in) :: has_concomitant !! True when concomitant coefficients should be stored.
      integer :: a, j, kk, component_df, concomitant_df
      call pack_common(result, active, posterior, logp, weights, first, coef_prior, has_concomitant)
      kk = result%k
      allocate(result%multinomial_coef(size(coef,1),size(coef,2),kk))
      a = 0
      do j = 1, size(active)
         if (.not. active(j)) cycle
         a = a + 1
         result%multinomial_coef(:,:,a) = coef(:,:,j)
      end do
      component_df = size(coef,1) * size(coef,2)
      if (has_concomitant) then
         concomitant_df = size(coef_prior,1) * max(0,kk - 1)
      else
         concomitant_df = max(0,kk - 1)
      end if
      result%df = kk * component_df + concomitant_df
   end subroutine pack_multinomial_result

   subroutine pack_multivariate_result(result, active, posterior, logp, center, covariance, probability, lambda_value, variance, &
                                       factor_loadings, uniqueness, marginal_variance, factors, weights, first, coef_prior, &
                                       has_concomitant, diagonal, truncated, binary_mask)
      type(flexmix_result), intent(inout) :: result !! Multivariate fit result receiving packed component parameters.
      logical, intent(in) :: active(:) !! Active-component mask over original components.
      real(dp), intent(in) :: posterior(:,:) !! Final posterior probabilities, shape `(n, k0)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized posterior log weights, shape `(n, k0)`.
      real(dp), intent(in) :: center(:,:) !! Component locations, shape `(d, k0)`.
      real(dp), intent(in) :: covariance(:,:,:) !! Component covariance matrices, shape `(d, d, k0)`.
      real(dp), intent(in) :: probability(:,:) !! Bernoulli probabilities, shape `(d, k0)`.
      real(dp), intent(in) :: lambda_value(:,:) !! Poisson means, shape `(d, k0)`.
      real(dp), intent(in) :: variance(:,:) !! Independent continuous variances for combined models, shape `(d, k0)`.
      real(dp), intent(in) :: factor_loadings(:,:,:) !! Correlation-scale factor loadings `(d,factors,k0)`; absent when unused.
      real(dp), intent(in) :: uniqueness(:,:) !! Correlation-scale factor uniquenesses, shape `(d,k0)`, or placeholders.
      real(dp), intent(in) :: marginal_variance(:,:) !! Unconstrained marginal variances, shape `(d,k0)`, or placeholders.
      integer, intent(in) :: factors !! Number of factor columns retained for factor-analyzer components.
      real(dp), intent(in) :: weights(:) !! Case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask, size `n`.
      real(dp), intent(in) :: coef_prior(:,:) !! Concomitant coefficient matrix or placeholder.
      logical, intent(in) :: has_concomitant !! True when concomitant coefficients should be stored.
      logical, intent(in) :: diagonal !! True when multivariate-normal covariance is constrained diagonal.
      logical, intent(in) :: truncated !! True when multivariate-binary likelihood is zero-truncated.
      logical, intent(in) :: binary_mask(:) !! Bernoulli-column mask used by the combined model.
      integer :: a, b, j, kk, m, component_df, concomitant_df, d
      call pack_common(result, active, posterior, logp, weights, first, coef_prior, has_concomitant)
      kk = result%k
      d = size(center,1)
      result%diagonal_covariance = diagonal
      result%truncated_binary = truncated
      select case (result%model_kind)
      case (flexmix_model_mvnorm)
         allocate(result%center(d,kk), result%covariance(d,d,kk))
         a = 0
         do j = 1, size(active)
            if (active(j)) then
               a = a + 1
               result%center(:,a) = center(:,j)
               result%covariance(:,:,a) = covariance(:,:,j)
            end if
         end do
         if (diagonal) then
            component_df = 2 * d
         else
            component_df = (3 * d + d * d) / 2
         end if
      case (flexmix_model_mvbinary)
         allocate(result%probability(d,kk))
         a = 0
         do j = 1, size(active)
            if (active(j)) then
               a = a + 1
               result%probability(:,a) = probability(:,j)
            end if
         end do
         component_df = d
      case (flexmix_model_mvpois)
         allocate(result%lambda(d,kk))
         a = 0
         do j = 1, size(active)
            if (active(j)) then
               a = a + 1
               result%lambda(:,a) = lambda_value(:,j)
            end if
         end do
         component_df = d
      case (flexmix_model_mvcombi)
         allocate(result%center(d,kk), result%covariance(d,d,kk), result%binary_mask(d))
         result%binary_mask = binary_mask
         result%covariance = 0.0_dp
         a = 0
         do j = 1, size(active)
            if (active(j)) then
               a = a + 1
               result%center(:,a) = center(:,j)
               result%covariance(:,:,a) = 0.0_dp
               m = 0
               do b = 1, d
                  if (.not. binary_mask(b)) then
                     m = m + 1
                     result%covariance(b,b,a) = variance(m,j)
                  end if
               end do
            end if
         end do
         component_df = d + count(.not. binary_mask)
      case (flexmix_model_factanal)
         result%factors = factors
         allocate(result%center(d,kk), result%covariance(d,d,kk), result%factor_loadings(d,factors,kk))
         allocate(result%uniqueness(d,kk), result%marginal_variance(d,kk))
         a = 0
         do j = 1, size(active)
            if (active(j)) then
               a = a + 1
               result%center(:,a) = center(:,j)
               result%covariance(:,:,a) = covariance(:,:,j)
               result%factor_loadings(:,:,a) = factor_loadings(:,:,j)
               result%uniqueness(:,a) = uniqueness(:,j)
               result%marginal_variance(:,a) = marginal_variance(:,j)
            end if
         end do
         component_df = (factors + 2) * d
      case default
         component_df = 0
      end select
      if (has_concomitant) then
         concomitant_df = size(coef_prior,1) * max(0,kk - 1)
      else
         concomitant_df = max(0,kk - 1)
      end if
      result%df = kk * component_df + concomitant_df
   end subroutine pack_multivariate_result

   subroutine pack_univariate_result(result, active, posterior, logp, par1, par2, weights, first, coef_prior, has_concomitant)
      type(flexmix_result), intent(inout) :: result !! Univariate fit result receiving packed distribution parameters.
      logical, intent(in) :: active(:) !! Active-component mask over original components.
      real(dp), intent(in) :: posterior(:,:) !! Final posterior probabilities, shape `(n, k0)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized posterior log weights, shape `(n, k0)`.
      real(dp), intent(in) :: par1(:) !! First fitted parameter for each original component.
      real(dp), intent(in) :: par2(:) !! Second fitted parameter for each original component when applicable.
      real(dp), intent(in) :: weights(:) !! Case weights, size `n`.
      logical, intent(in) :: first(:) !! Group-first mask, size `n`.
      real(dp), intent(in) :: coef_prior(:,:) !! Concomitant coefficient matrix or placeholder.
      logical, intent(in) :: has_concomitant !! True when concomitant coefficients should be stored.
      integer :: a, j, kk, component_df, concomitant_df
      call pack_common(result, active, posterior, logp, weights, first, coef_prior, has_concomitant)
      kk = result%k
      select case (result%model_kind)
      case (flexmix_model_lognormal)
         allocate(result%center(1,kk), result%sigma(kk))
         component_df = 2
      case (flexmix_model_exponential)
         allocate(result%rate(kk))
         component_df = 1
      case (flexmix_model_inverse_gaussian)
         allocate(result%center(1,kk), result%lambda(1,kk))
         component_df = 2
      case (flexmix_model_gamma)
         allocate(result%shape(kk), result%rate(kk))
         component_df = 2
      case (flexmix_model_weibull)
         allocate(result%shape(kk), result%scale(kk))
         component_df = 2
      case default
         component_df = 0
      end select
      a = 0
      do j = 1, size(active)
         if (.not. active(j)) cycle
         a = a + 1
         select case (result%model_kind)
         case (flexmix_model_lognormal)
            result%center(1,a) = par1(j)
            result%sigma(a) = par2(j)
         case (flexmix_model_exponential)
            result%rate(a) = par1(j)
         case (flexmix_model_inverse_gaussian)
            result%center(1,a) = par1(j)
            result%lambda(1,a) = par2(j)
         case (flexmix_model_gamma)
            result%shape(a) = par1(j)
            result%rate(a) = par2(j)
         case (flexmix_model_weibull)
            result%shape(a) = par1(j)
            result%scale(a) = par2(j)
         end select
      end do
      if (has_concomitant) then
         concomitant_df = size(coef_prior,1) * max(0,kk - 1)
      else
         concomitant_df = max(0,kk - 1)
      end if
      result%df = kk * component_df + concomitant_df
   end subroutine pack_univariate_result


end module flexmix_core
