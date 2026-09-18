! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_api
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_control, flexmix_result, flexmix_step_result
   use flexmix_types, only : flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial
   use flexmix_types, only : flexmix_model_mvnorm, flexmix_model_mvbinary, flexmix_model_mvpois, flexmix_model_mvcombi
   use flexmix_types, only : flexmix_model_lognormal, flexmix_model_exponential, flexmix_model_inverse_gaussian
   use flexmix_types, only : flexmix_model_gamma, flexmix_model_weibull
   use flexmix_types, only : flexmix_model_gamma_regression, flexmix_model_multinomial
   use flexmix_types, only : flexmix_model_ziglm_poisson, flexmix_model_ziglm_binomial
   use flexmix_types, only : flexmix_model_robust_gaussian, flexmix_model_robust_poisson, flexmix_model_factanal
   use flexmix_types, only : flexmix_model_conditional_logit, flexmix_model_lmm, flexmix_model_lmer, flexmix_model_lmmc, &
                             flexmix_model_lmc
   use flexmix_types, only : flexmix_role_regular, flexmix_role_structural_zero, flexmix_role_robust_background
   use flexmix_core, only : flexmix_fit_regression, flexmix_fit_binomial, flexmix_fit_multinomial, flexmix_fit_glmnet
   use flexmix_core, only : flexmix_fit_mgcv
   use flexmix_core, only : flexmix_fit_conditional_logit
   use flexmix_core, only : flexmix_fit_ziglm_poisson, flexmix_fit_ziglm_binomial
   use flexmix_core, only : flexmix_fit_robust_gaussian, flexmix_fit_robust_poisson
   use flexmix_core, only : flexmix_fit_multivariate, flexmix_fit_univariate
   use flexmix_numeric, only : inverse_logdet_spd, clip_probability
   use flexmix_components, only : multinomial_regression_probabilities
   implicit none
   private
   public :: flexmix_gaussian
   public :: flexmix_poisson
   public :: flexmix_binomial
   public :: flexmix_gamma
   public :: flexmix_glmnet_gaussian
   public :: flexmix_glmnet_poisson
   public :: flexmix_glmnet_binomial
   public :: flexmix_mgcv_gaussian
   public :: flexmix_mgcv_poisson
   public :: flexmix_mgcv_binomial
   public :: flexmix_multinomial
   public :: flexmix_conditional_logit
   public :: flexmix_ziglm_poisson
   public :: flexmix_ziglm_binomial
   public :: flexmix_robust_gaussian
   public :: flexmix_robust_poisson
   public :: flexmix_mvnorm
   public :: flexmix_factanal
   public :: flexmix_mvbinary
   public :: flexmix_mvpois
   public :: flexmix_mvcombi
   public :: flexmix_dist1
   public :: flexmix_norm1
   public :: flexmix_loglik
   public :: flexmix_nobs
   public :: flexmix_aic
   public :: flexmix_bic
   public :: flexmix_cloglik
   public :: flexmix_icl
   public :: flexmix_eic
   public :: flexmix_get_k
   public :: flexmix_get_obs
   public :: flexmix_get_prior
   public :: flexmix_get_posterior
   public :: flexmix_get_clusters
   public :: flexmix_get_parameters
   public :: flexmix_predict_regression
   public :: flexmix_predict_multinomial
   public :: flexmix_relabel
   public :: kl_divergence_matrix
   public :: kl_divergence_regression
   public :: kl_divergence_mvnorm
   public :: step_flexmix_gaussian
   public :: step_flexmix_unique
   public :: get_model
   public :: flexmix_check_result
   public :: flexmix_remove_component
   public :: flexmix_refit_gaussian
   public :: flexmix_simulate_regression
   public :: flexmix_make_distribution
   public :: ex_linear
   public :: ex_npreg
   public :: ex_nclus

contains

   subroutine flexmix_gaussian(x, y, k, result, control, case_weights, group, concomitant_x, &
                               initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Numeric Gaussian-regression design matrix, shape `(n, p)`, including any intercept.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted Gaussian mixture-of-regressions result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive identity-link offset, size `n`.
      call flexmix_fit_regression(x, y, k, flexmix_model_gaussian, result, control, case_weights, group, &
                                 concomitant_x, initial_cluster, initial_posterior, offset)
   end subroutine flexmix_gaussian

   subroutine flexmix_poisson(x, y, k, result, control, case_weights, group, concomitant_x, &
                              initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Numeric Poisson-regression design matrix, shape `(n, p)`, including any intercept.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson count response, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted Poisson mixture-of-regressions result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive log-link offset, size `n`.
      call flexmix_fit_regression(x, y, k, flexmix_model_poisson, result, control, case_weights, group, &
                                 concomitant_x, initial_cluster, initial_posterior, offset)
   end subroutine flexmix_poisson

   subroutine flexmix_binomial(x, success, failure, k, result, control, case_weights, group, concomitant_x, &
                               initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Numeric binomial-regression design matrix, shape `(n, p)`, including any intercept.
      real(dp), intent(in) :: success(:) !! Nonnegative success counts, size `n`.
      real(dp), intent(in) :: failure(:) !! Nonnegative failure counts, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted binomial mixture-of-regressions result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive logit-link offset, size `n`.
      call flexmix_fit_binomial(x, success, failure, k, result, control, case_weights, group, concomitant_x, &
                               initial_cluster, initial_posterior, offset)
   end subroutine flexmix_binomial

   subroutine flexmix_glmnet_gaussian(x, y, k, result, control, case_weights, group, concomitant_x, initial_cluster, &
                                           initial_posterior, offset, adaptive, select, alpha, nfolds, fold_id, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Gaussian penalized-regression design matrix `(n,p)` with an intercept in column one.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted adaptive-lasso/elastic-net Gaussian mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive identity-link offset, size `n`.
      logical, intent(in), optional :: adaptive !! Use adaptive penalty factors from weighted unpenalized fits; default true.
      logical, intent(in), optional :: select(:) !! Penalized-term mask for non-intercept columns, size `p-1`.
      real(dp), intent(in), optional :: alpha !! Elastic-net mixing fraction in `[0,1]`; default one.
      integer, intent(in), optional :: nfolds !! Number of deterministic CV folds when `fold_id` is absent; default five.
      integer, intent(in), optional :: fold_id(:) !! Optional fixed one-based CV fold labels, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative penalty grid searched by cross-validation.
      call flexmix_fit_glmnet(x, y, k, flexmix_model_gaussian, result, control=control, case_weights=case_weights, &
                              group=group, concomitant_x=concomitant_x, initial_cluster=initial_cluster, &
                              initial_posterior=initial_posterior, offset=offset, adaptive=adaptive, select=select, &
                              alpha=alpha, nfolds=nfolds, fold_id=fold_id, lambda_grid=lambda_grid)
   end subroutine flexmix_glmnet_gaussian

   subroutine flexmix_glmnet_poisson(x, y, k, result, control, case_weights, group, concomitant_x, initial_cluster, &
                                    initial_posterior, offset, adaptive, select, alpha, nfolds, fold_id, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Poisson penalized-regression design matrix `(n,p)` with an intercept in column one.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson response vector, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted adaptive-lasso/elastic-net Poisson mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive log-link offset, size `n`.
      logical, intent(in), optional :: adaptive !! Use adaptive penalty factors from weighted unpenalized fits; default true.
      logical, intent(in), optional :: select(:) !! Penalized-term mask for non-intercept columns, size `p-1`.
      real(dp), intent(in), optional :: alpha !! Elastic-net mixing fraction in `[0,1]`; default one.
      integer, intent(in), optional :: nfolds !! Number of deterministic CV folds when `fold_id` is absent; default five.
      integer, intent(in), optional :: fold_id(:) !! Optional fixed one-based CV fold labels, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative penalty grid searched by cross-validation.
      call flexmix_fit_glmnet(x, y, k, flexmix_model_poisson, result, control=control, case_weights=case_weights, &
                              group=group, concomitant_x=concomitant_x, initial_cluster=initial_cluster, &
                              initial_posterior=initial_posterior, offset=offset, adaptive=adaptive, select=select, &
                              alpha=alpha, nfolds=nfolds, fold_id=fold_id, lambda_grid=lambda_grid)
   end subroutine flexmix_glmnet_poisson

   subroutine flexmix_glmnet_binomial(x, success, failure, k, result, control, case_weights, group, concomitant_x, &
                                     initial_cluster, initial_posterior, offset, adaptive, select, alpha, nfolds, &
                                     fold_id, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Binomial penalized-regression design matrix `(n,p)` with an intercept in column one.
      real(dp), intent(in) :: success(:) !! Nonnegative binomial success counts, size `n`.
      real(dp), intent(in) :: failure(:) !! Nonnegative binomial failure counts, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted adaptive-lasso/elastic-net binomial mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional group labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive logit offset, size `n`.
      logical, intent(in), optional :: adaptive !! Use adaptive penalty factors from weighted unpenalized fits; default true.
      logical, intent(in), optional :: select(:) !! Penalized-term mask for non-intercept columns, size `p-1`.
      real(dp), intent(in), optional :: alpha !! Elastic-net mixing fraction in `[0,1]`; default one.
      integer, intent(in), optional :: nfolds !! Number of deterministic CV folds when `fold_id` is absent; default five.
      integer, intent(in), optional :: fold_id(:) !! Optional fixed one-based CV fold labels, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative penalty grid searched by cross-validation.
      real(dp), allocatable :: trials(:)
      allocate(trials(size(success)))
      trials = success + failure
      call flexmix_fit_glmnet(x, success, k, flexmix_model_binomial, result, trials=trials, control=control, &
                              case_weights=case_weights, group=group, concomitant_x=concomitant_x, &
                              initial_cluster=initial_cluster, initial_posterior=initial_posterior, offset=offset, &
                              adaptive=adaptive, select=select, alpha=alpha, nfolds=nfolds, fold_id=fold_id, &
                              lambda_grid=lambda_grid)
   end subroutine flexmix_glmnet_binomial

   subroutine flexmix_mgcv_gaussian(x, y, penalty_matrix, k, result, control, case_weights, group, concomitant_x, &
                                        initial_cluster, initial_posterior, offset, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Prefit Gaussian GAM/design matrix `(n,p)` containing parametric and smooth-basis.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: penalty_matrix(:,:) !! Symmetric smoothing penalty matrix `(p,p)`; zero rows/columns are unpenalized.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted Gaussian smooth-regression mixture with smoothing diagnostics.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive identity-link offset, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative smoothing grid; generated when absent.
      call flexmix_fit_mgcv(x, y, k, flexmix_model_gaussian, penalty_matrix, result, control=control, &
                            case_weights=case_weights, group=group, concomitant_x=concomitant_x, &
                            initial_cluster=initial_cluster, initial_posterior=initial_posterior, offset=offset, &
                            lambda_grid=lambda_grid)
   end subroutine flexmix_mgcv_gaussian

   subroutine flexmix_mgcv_poisson(x, y, penalty_matrix, k, result, control, case_weights, group, concomitant_x, &
                                  initial_cluster, initial_posterior, offset, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Prefit Poisson GAM/design matrix `(n,p)` containing parametric and smooth-basis columns.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson response vector, size `n`.
      real(dp), intent(in) :: penalty_matrix(:,:) !! Symmetric smoothing penalty matrix `(p,p)`; zero rows/columns are unpenalized.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted Poisson smooth-regression mixture with smoothing diagnostics.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive log-link offset, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative smoothing grid; generated when absent.
      call flexmix_fit_mgcv(x, y, k, flexmix_model_poisson, penalty_matrix, result, control=control, &
                            case_weights=case_weights, group=group, concomitant_x=concomitant_x, &
                            initial_cluster=initial_cluster, initial_posterior=initial_posterior, offset=offset, &
                            lambda_grid=lambda_grid)
   end subroutine flexmix_mgcv_poisson

   subroutine flexmix_mgcv_binomial(x, success, failure, penalty_matrix, k, result, control, case_weights, group, &
                                   concomitant_x, initial_cluster, initial_posterior, offset, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Prefit binomial GAM/design matrix `(n,p)` containing parametric and smooth-basis.
      real(dp), intent(in) :: success(:) !! Nonnegative binomial success counts, size `n`.
      real(dp), intent(in) :: failure(:) !! Nonnegative binomial failure counts, size `n`.
      real(dp), intent(in) :: penalty_matrix(:,:) !! Symmetric smoothing penalty matrix `(p,p)`; zero rows/columns are unpenalized.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted binomial smooth-regression mixture with smoothing diagnostics.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive logit offset, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative smoothing grid; generated when absent.
      real(dp), allocatable :: trials(:)
      allocate(trials(size(success)))
      trials = success + failure
      call flexmix_fit_mgcv(x, success, k, flexmix_model_binomial, penalty_matrix, result, trials=trials, &
                            control=control, case_weights=case_weights, group=group, concomitant_x=concomitant_x, &
                            initial_cluster=initial_cluster, initial_posterior=initial_posterior, offset=offset, &
                            lambda_grid=lambda_grid)
   end subroutine flexmix_mgcv_binomial

   subroutine flexmix_gamma(x, y, k, result, control, case_weights, group, concomitant_x, &
                            initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Numeric Gamma-regression design matrix, shape `(n, p)`, including any intercept column.
      real(dp), intent(in) :: y(:) !! Strictly positive Gamma response vector, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted Gamma mixture-of-regressions result using R's default inverse link.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive inverse-link offset, size `n`.
      call flexmix_fit_regression(x, y, k, flexmix_model_gamma_regression, result, control, case_weights, group, &
                                 concomitant_x, initial_cluster, initial_posterior, offset)
   end subroutine flexmix_gamma

   subroutine flexmix_multinomial(x, y_class, k, result, control, case_weights, group, concomitant_x, &
                                  initial_cluster, initial_posterior)
      real(dp), intent(in) :: x(:,:) !! Multinomial-regression design matrix, shape `(n, p)`, including an intercept if wanted.
      integer, intent(in) :: y_class(:) !! One-based response class labels covering all integers from one through the maximum.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted mixture of multinomial logistic regressions.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial mixture-component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial mixture posterior matrix, shape `(n, k)`.
      call flexmix_fit_multinomial(x, y_class, k, result, control, case_weights, group, concomitant_x, &
                                  initial_cluster, initial_posterior)
   end subroutine flexmix_multinomial

   subroutine flexmix_conditional_logit(x, y, strata, k, result, control, case_weights, &
                                          initial_cluster, initial_posterior)
      real(dp), intent(in) :: x(:,:) !! Conditional-logit covariate matrix without an intercept, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Binary selected-alternative indicator, size `n`, with one selected row per stratum.
      integer, intent(in) :: strata(:) !! Integer choice-set labels, size `n`; rows sharing a label form one stratum.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted mixture of conditional-logit components.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative stratum-frequency weights, constant within.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting mixture labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior matrix, shape `(n,k)`.
      call flexmix_fit_conditional_logit(x, y, strata, k, result, control, case_weights, &
                                        initial_cluster, initial_posterior)
   end subroutine flexmix_conditional_logit

   subroutine flexmix_ziglm_poisson(x, y, k, result, control, case_weights, group, concomitant_x, &
                                      initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Poisson design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson counts, size `n`; exact zeros are eligible for the structural.
      integer, intent(in) :: k !! Requested component count including the structural-zero component.
      type(flexmix_result), intent(out) :: result !! Fitted zero-inflated Poisson FlexMix-style result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, classification, and minprior settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive Poisson log-link offset, size `n`.
      call flexmix_fit_ziglm_poisson(x, y, k, result, control, case_weights, group, concomitant_x, &
                                     initial_cluster, initial_posterior, offset)
   end subroutine flexmix_ziglm_poisson

   subroutine flexmix_ziglm_binomial(x, success, failure, k, result, control, case_weights, group, &
                                     concomitant_x, initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Binomial design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: success(:) !! Nonnegative success counts, size `n`; zero-success rows fit the structural.
      real(dp), intent(in) :: failure(:) !! Nonnegative failure counts, size `n`, paired with `success`.
      integer, intent(in) :: k !! Requested component count including the structural-zero component.
      type(flexmix_result), intent(out) :: result !! Fitted zero-inflated binomial FlexMix-style result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, classification, and minprior settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive binomial logit-link offset, size `n`.
      call flexmix_fit_ziglm_binomial(x, success, failure, k, result, control, case_weights, group, concomitant_x, &
                                      initial_cluster, initial_posterior, offset)
   end subroutine flexmix_ziglm_binomial

   subroutine flexmix_robust_gaussian(x, y, k, result, bgw, control, case_weights, group, concomitant_x, &
                                      initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Gaussian design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      integer, intent(in) :: k !! Requested component count including the robust background component.
      type(flexmix_result), intent(out) :: result !! Fitted robust Gaussian FlexMix-style result.
      logical, intent(in), optional :: bgw !! Use current posterior weights for the background component when true; default.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, classification, and minprior settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive Gaussian identity-link offset, size `n`.
      call flexmix_fit_robust_gaussian(x, y, k, result, bgw, control, case_weights, group, concomitant_x, &
                                       initial_cluster, initial_posterior, offset)
   end subroutine flexmix_robust_gaussian

   subroutine flexmix_robust_poisson(x, y, k, result, bgw, control, case_weights, group, concomitant_x, &
                                     initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:) !! Poisson design matrix with an intercept in column one, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson count response, size `n`.
      integer, intent(in) :: k !! Requested component count including the robust background component.
      type(flexmix_result), intent(out) :: result !! Fitted robust Poisson FlexMix-style result.
      logical, intent(in), optional :: bgw !! Use current posterior weights for the background component when true; default.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, classification, and minprior settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior probabilities, shape `(n, k)`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive Poisson log-link offset, size `n`.
      call flexmix_fit_robust_poisson(x, y, k, result, bgw, control, case_weights, group, concomitant_x, &
                                      initial_cluster, initial_posterior, offset)
   end subroutine flexmix_robust_poisson

   subroutine flexmix_mvnorm(y, k, result, diagonal, control, case_weights, group, concomitant_x, &
                             initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:,:) !! Numeric clustering matrix, shape `(n, d)`.
      integer, intent(in) :: k !! Requested initial number of Gaussian components.
      type(flexmix_result), intent(out) :: result !! Fitted multivariate-normal mixture result.
      logical, intent(in), optional :: diagonal !! Use diagonal covariance matrices when true; default is full covariance.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      call flexmix_fit_multivariate(y, k, flexmix_model_mvnorm, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior, diagonal=diagonal)
   end subroutine flexmix_mvnorm

   subroutine flexmix_factanal(y, k, factors, result, control, case_weights, group, concomitant_x, &
                                initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:,:) !! Multivariate numeric response matrix, shape `(n,d)`, modeled by factor-analyzer.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      integer, intent(in) :: factors !! Number of latent factors per component, subject to the `factanal` degrees-of-freedom.
      type(flexmix_result), intent(out) :: result !! Fitted mixture-of-factor-analyzers result with centers, covariances.
      type(flexmix_control), intent(in), optional :: control !! Optional EM convergence, minprior, and classification settings.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional grouped-observation labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional multinomial-logit concomitant design matrix `(n,q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior matrix, shape `(n,k)`.
      call flexmix_fit_multivariate(y, k, flexmix_model_factanal, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior, factors=factors)
   end subroutine flexmix_factanal

   subroutine flexmix_norm1(y, k, result, control, case_weights, group, concomitant_x, initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:) !! Univariate normal observations, size `n`.
      integer, intent(in) :: k !! Requested initial number of normal components.
      type(flexmix_result), intent(out) :: result !! Fitted one-dimensional normal mixture stored as normal components.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      real(dp), allocatable :: yy(:,:)
      allocate(yy(size(y),1))
      yy(:,1) = y
      call flexmix_fit_multivariate(yy, k, flexmix_model_mvnorm, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior, diagonal=.true.)
   end subroutine flexmix_norm1

   subroutine flexmix_mvbinary(y, k, result, truncated, control, case_weights, group, concomitant_x, &
                               initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:,:) !! Binary response matrix containing zeros and ones, shape `(n, d)`.
      integer, intent(in) :: k !! Requested initial number of Bernoulli-product components.
      type(flexmix_result), intent(out) :: result !! Fitted multivariate-binary mixture result.
      logical, intent(in), optional :: truncated !! Condition the Bernoulli product on at least one success when true.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      call flexmix_fit_multivariate(y, k, flexmix_model_mvbinary, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior, truncated=truncated)
   end subroutine flexmix_mvbinary

   subroutine flexmix_mvpois(y, k, result, control, case_weights, group, concomitant_x, initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:,:) !! Nonnegative multivariate count matrix, shape `(n, d)`.
      integer, intent(in) :: k !! Requested initial number of independent-Poisson components.
      type(flexmix_result), intent(out) :: result !! Fitted multivariate-Poisson mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      call flexmix_fit_multivariate(y, k, flexmix_model_mvpois, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior)
   end subroutine flexmix_mvpois

   subroutine flexmix_mvcombi(y, binary, k, result, control, case_weights, group, concomitant_x, &
                              initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:,:) !! Mixed binary/continuous response matrix, shape `(n, d)`.
      logical, intent(in) :: binary(:) !! Mask of Bernoulli columns, size `d`; false columns use independent normal densities.
      integer, intent(in) :: k !! Requested initial number of combined components.
      type(flexmix_result), intent(out) :: result !! Fitted combined binary/continuous mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      call flexmix_fit_multivariate(y, k, flexmix_model_mvcombi, result, control, case_weights, group, concomitant_x, &
                                   initial_cluster, initial_posterior, binary=binary)
   end subroutine flexmix_mvcombi

   subroutine flexmix_dist1(y, k, distribution, result, control, case_weights, group, concomitant_x, &
                            initial_cluster, initial_posterior)
      real(dp), intent(in) :: y(:) !! Positive univariate observations, size `n`.
      integer, intent(in) :: k !! Requested initial number of distribution components.
      character(len=*), intent(in) :: distribution !! Distribution name: `lnorm`, `exp`, `invgauss`, `gamma`, or `weibull`.
      type(flexmix_result), intent(out) :: result !! Fitted univariate distribution mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional FlexMix-style EM controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional integer grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant softmax design matrix, shape `(n, q)`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial posterior matrix, shape `(n, k)`.
      integer :: model_kind
      select case (trim(adjustl(distribution)))
      case ('lnorm', 'lognormal')
         model_kind = flexmix_model_lognormal
      case ('exp', 'exponential')
         model_kind = flexmix_model_exponential
      case ('invgauss', 'inverse_gaussian')
         model_kind = flexmix_model_inverse_gaussian
      case ('gamma')
         model_kind = flexmix_model_gamma
      case ('weibull')
         model_kind = flexmix_model_weibull
      case default
         result%status = -20
         result%model_kind = 0
         return
      end select
      call flexmix_fit_univariate(y, k, model_kind, result, control, case_weights, group, concomitant_x, &
                                  initial_cluster, initial_posterior)
   end subroutine flexmix_dist1

   pure function flexmix_loglik(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result whose stored maximized log likelihood is returned.
      real(dp) :: value
      value = result%loglik
   end function flexmix_loglik

   pure function flexmix_nobs(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result containing original case weights when available.
      real(dp) :: value
      if (allocated(result%case_weights)) then
         value = sum(result%case_weights)
      else
         value = real(result%n, dp)
      end if
   end function flexmix_nobs

   pure function flexmix_aic(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture supplying log likelihood and effective parameter count.
      real(dp) :: value
      if (result%effective_df > 0.0_dp) then
         value = -2.0_dp * result%loglik + 2.0_dp * result%effective_df
      else
         value = -2.0_dp * result%loglik + 2.0_dp * real(result%df, dp)
      end if
   end function flexmix_aic

   pure function flexmix_bic(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture supplying log likelihood, parameter count, and sample size.
      real(dp) :: value
      if (result%effective_df > 0.0_dp) then
         value = -2.0_dp * result%loglik + result%effective_df * log(max(1.0_dp, flexmix_nobs(result)))
      else
         value = -2.0_dp * result%loglik + real(result%df, dp) * log(max(1.0_dp, flexmix_nobs(result)))
      end if
   end function flexmix_bic

   pure function flexmix_cloglik(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture with unscaled posterior logs and modal labels.
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      if (.not. allocated(result%log_posterior_unscaled)) return
      do i = 1, result%n
         if (allocated(result%group_first)) then
            if (.not. result%group_first(i)) cycle
         end if
         value = value + result%log_posterior_unscaled(i, result%cluster(i))
      end do
   end function flexmix_cloglik

   pure function flexmix_icl(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture for the classification-likelihood information criterion.
      real(dp) :: value
      if (result%effective_df > 0.0_dp) then
         value = -2.0_dp * flexmix_cloglik(result) + result%effective_df * log(max(1.0_dp, flexmix_nobs(result)))
      else
         value = -2.0_dp * flexmix_cloglik(result) + real(result%df, dp) * log(max(1.0_dp, flexmix_nobs(result)))
      end if
   end function flexmix_icl

   pure function flexmix_eic(result) result(value)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result whose group-first posterior entropy is summarized.
      real(dp) :: value
      real(dp) :: p, total
      integer :: i, j, nfirst
      if (result%k <= 1 .or. .not. allocated(result%posterior)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      total = 0.0_dp
      nfirst = 0
      do i = 1, result%n
         if (allocated(result%group_first)) then
            if (.not. result%group_first(i)) cycle
         end if
         nfirst = nfirst + 1
         do j = 1, result%k
            p = result%posterior(i,j)
            if (p > 0.0_dp) total = total + p * max(-1000.0_dp, log(p))
         end do
      end do
      if (nfirst == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = 1.0_dp + total / (real(nfirst,dp) * log(real(result%k,dp)))
      end if
   end function flexmix_eic

   pure function flexmix_get_k(result) result(k)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result whose retained component count is requested.
      integer :: k
      k = result%k
   end function flexmix_get_k

   pure function flexmix_get_obs(result) result(n)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result whose stored row count is requested.
      integer :: n
      n = result%n
   end function flexmix_get_obs

   subroutine flexmix_get_prior(result, prior)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result containing marginal component priors.
      real(dp), allocatable, intent(out) :: prior(:) !! Copy of the retained marginal prior vector.
      if (allocated(result%prior)) then
         allocate(prior(size(result%prior)))
         prior = result%prior
      else
         allocate(prior(0))
      end if
   end subroutine flexmix_get_prior

   subroutine flexmix_get_posterior(result, posterior, unscaled)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result containing scaled and unscaled posterior matrices.
      real(dp), allocatable, intent(out) :: posterior(:,:) !! Copy of the requested posterior matrix.
      logical, intent(in), optional :: unscaled !! Return unnormalized posterior likelihood weights when true; default false.
      logical :: raw
      raw = .false.
      if (present(unscaled)) raw = unscaled
      if (raw .and. allocated(result%posterior_unscaled)) then
         allocate(posterior(size(result%posterior_unscaled,1), size(result%posterior_unscaled,2)))
         posterior = result%posterior_unscaled
      else if (allocated(result%posterior)) then
         allocate(posterior(size(result%posterior,1), size(result%posterior,2)))
         posterior = result%posterior
      else
         allocate(posterior(0,0))
      end if
   end subroutine flexmix_get_posterior

   subroutine flexmix_get_clusters(result, cluster)
      type(flexmix_result), intent(in) :: result !! Fitted mixture result containing modal component assignments.
      integer, allocatable, intent(out) :: cluster(:) !! Copy of one-based modal component labels, size `n`.
      if (allocated(result%cluster)) then
         allocate(cluster(size(result%cluster)))
         cluster = result%cluster
      else
         allocate(cluster(0))
      end if
   end subroutine flexmix_get_clusters

   subroutine flexmix_get_parameters(result, parameters, info)
      type(flexmix_result), intent(in) :: result !! Fitted mixture whose component parameters are requested numerically.
      real(dp), allocatable, intent(out) :: parameters(:,:) !! Columns are components; rows contain family-specific parameters.
      integer, intent(out) :: info !! Zero on success; nonzero when the model family has no generic numeric packing here.
      integer :: a, d, p, rows
      info = 0
      select case (result%model_kind)
      case (flexmix_model_gaussian, flexmix_model_robust_gaussian, flexmix_model_lmc)
         p = size(result%beta,1)
         allocate(parameters(p+1,result%k))
         parameters(1:p,:) = result%beta
         parameters(p+1,:) = result%sigma
      case (flexmix_model_poisson, flexmix_model_binomial, flexmix_model_ziglm_poisson, &
            flexmix_model_ziglm_binomial, flexmix_model_robust_poisson, flexmix_model_conditional_logit)
         allocate(parameters(size(result%beta,1),result%k))
         parameters = result%beta
      case (flexmix_model_gamma_regression)
         p = size(result%beta,1)
         allocate(parameters(p+1,result%k))
         parameters(1:p,:) = result%beta
         parameters(p+1,:) = result%shape
      case (flexmix_model_lmm, flexmix_model_lmer, flexmix_model_lmmc)
         p = size(result%beta,1)
         d = size(result%random_covariance,1)
         rows = p + d*d + 1
         allocate(parameters(rows,result%k))
         do a = 1, result%k
            parameters(1:p,a) = result%beta(:,a)
            parameters(p+1:p+d*d,a) = reshape(result%random_covariance(:,:,a), [d*d])
            parameters(rows,a) = result%residual_variance(a)
         end do
      case (flexmix_model_multinomial)
         rows = size(result%multinomial_coef,1) * size(result%multinomial_coef,2)
         allocate(parameters(rows,result%k))
         do a = 1, result%k
            parameters(:,a) = reshape(result%multinomial_coef(:,:,a), [rows])
         end do
      case (flexmix_model_mvnorm)
         d = result%d
         rows = d + d * d
         allocate(parameters(rows,result%k))
         do a = 1, result%k
            parameters(1:d,a) = result%center(:,a)
            parameters(d+1:rows,a) = reshape(result%covariance(:,:,a), [d*d])
         end do
      case (flexmix_model_factanal)
         d = result%d
         rows = 3*d + d*result%factors
         allocate(parameters(rows,result%k))
         do a = 1, result%k
            parameters(1:d,a) = result%center(:,a)
            parameters(d+1:2*d,a) = result%marginal_variance(:,a)
            parameters(2*d+1:2*d+d*result%factors,a) = reshape(result%factor_loadings(:,:,a), [d*result%factors])
            parameters(2*d+d*result%factors+1:rows,a) = result%uniqueness(:,a)
         end do
      case (flexmix_model_mvbinary)
         allocate(parameters(result%d,result%k))
         parameters = result%probability
      case (flexmix_model_mvpois)
         allocate(parameters(result%d,result%k))
         parameters = result%lambda
      case (flexmix_model_mvcombi)
         d = result%d
         allocate(parameters(2*d,result%k))
         parameters(1:d,:) = result%center
         do a = 1, result%k
            do p = 1, d
               parameters(d+p,a) = result%covariance(p,p,a)
            end do
         end do
      case (flexmix_model_lognormal)
         allocate(parameters(2,result%k))
         parameters(1,:) = result%center(1,:)
         parameters(2,:) = result%sigma
      case (flexmix_model_exponential)
         allocate(parameters(1,result%k))
         parameters(1,:) = result%rate
      case (flexmix_model_inverse_gaussian)
         allocate(parameters(2,result%k))
         parameters(1,:) = result%center(1,:)
         parameters(2,:) = result%lambda(1,:)
      case (flexmix_model_gamma)
         allocate(parameters(2,result%k))
         parameters(1,:) = result%shape
         parameters(2,:) = result%rate
      case (flexmix_model_weibull)
         allocate(parameters(2,result%k))
         parameters(1,:) = result%shape
         parameters(2,:) = result%scale
      case default
         allocate(parameters(0,0))
         info = 1
      end select
   end subroutine flexmix_get_parameters

   subroutine flexmix_predict_regression(result, x, prediction, info, offset)
      type(flexmix_result), intent(in) :: result !! Fitted Gaussian, Poisson, binomial, or Gamma mixture-of-regressions result.
      real(dp), intent(in) :: x(:,:) !! New numeric design matrix, shape `(nnew, p)` matching the fitted coefficient count.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Component-wise conditional means/probabilities, shape `(nnew.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported model kind or incompatible design/offset shape.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive family-scale offset for the new rows, size `nnew`.
      real(dp), allocatable :: eta(:)
      integer :: j
      if (.not. allocated(result%beta) .or. size(x,2) /= result%p) then
         allocate(prediction(0,0))
         info = 1
         return
      end if
      if (present(offset)) then
         if (size(offset) /= size(x,1)) then
            allocate(prediction(0,0))
            info = 3
            return
         end if
      end if
      allocate(prediction(size(x,1),result%k), eta(size(x,1)))
      do j = 1, result%k
         eta = matmul(x, result%beta(:,j))
         if (present(offset)) eta = eta + offset
         select case (result%model_kind)
         case (flexmix_model_gaussian, flexmix_model_robust_gaussian, flexmix_model_lmc, flexmix_model_lmm, &
               flexmix_model_lmer, flexmix_model_lmmc)
            prediction(:,j) = eta
         case (flexmix_model_poisson, flexmix_model_robust_poisson)
            prediction(:,j) = exp(max(-30.0_dp, min(30.0_dp, eta)))
         case (flexmix_model_ziglm_poisson)
            prediction(:,j) = exp(max(-30.0_dp, min(30.0_dp, eta)))
            if (allocated(result%component_role)) then
               if (result%component_role(j) == flexmix_role_structural_zero) prediction(:,j) = 0.0_dp
            end if
         case (flexmix_model_binomial)
            prediction(:,j) = 1.0_dp / (1.0_dp + exp(-max(-30.0_dp, min(30.0_dp, eta))))
         case (flexmix_model_ziglm_binomial)
            prediction(:,j) = 1.0_dp / (1.0_dp + exp(-max(-30.0_dp, min(30.0_dp, eta))))
            if (allocated(result%component_role)) then
               if (result%component_role(j) == flexmix_role_structural_zero) prediction(:,j) = 0.0_dp
            end if
         case (flexmix_model_gamma_regression)
            prediction(:,j) = 1.0_dp / max(1.0e-8_dp, min(1.0e8_dp, eta))
         case (flexmix_model_conditional_logit)
            prediction(:,j) = eta
         case default
            prediction = 0.0_dp
            info = 2
            return
         end select
      end do
      info = 0
   end subroutine flexmix_predict_regression

   subroutine flexmix_predict_multinomial(result, x, prediction, info)
      type(flexmix_result), intent(in) :: result !! Fitted multinomial-regression mixture result.
      real(dp), intent(in) :: x(:,:) !! New numeric design matrix, shape `(nnew, p)` matching the fitted model.
      real(dp), allocatable, intent(out) :: prediction(:,:,:) !! Class probabilities, shape `(nnew, nclass, k)`.
      integer, intent(out) :: info !! Zero on success; nonzero for the wrong model family or incompatible design shape.
      integer :: j
      if (result%model_kind /= flexmix_model_multinomial .or. .not. allocated(result%multinomial_coef) .or. &
          size(x,2) /= result%p) then
         allocate(prediction(0,0,0))
         info = 1
         return
      end if
      allocate(prediction(size(x,1),result%nclass,result%k))
      do j = 1, result%k
         call multinomial_regression_probabilities(x, result%multinomial_coef(:,:,j), prediction(:,:,j))
      end do
      info = 0
   end subroutine flexmix_predict_multinomial

   subroutine flexmix_relabel(result, permutation, info)
      type(flexmix_result), intent(inout) :: result !! Fitted result whose component order and labels are permuted in place.
      integer, intent(in) :: permutation(:) !! New-order vector of old component indices; each component must appear once.
      integer, intent(out) :: info !! Zero on success; nonzero when `permutation` is invalid.
      type(flexmix_result) :: old
      integer, allocatable :: inverse(:)
      integer :: i
      if (size(permutation) /= result%k) then
         info = 1
         return
      end if
      allocate(inverse(result%k))
      inverse = 0
      do i = 1, result%k
         if (permutation(i) < 1 .or. permutation(i) > result%k) then
            info = 2
            return
         end if
         if (inverse(permutation(i)) /= 0) then
            info = 3
            return
         end if
         inverse(permutation(i)) = i
      end do
      old = result
      if (allocated(old%prior)) result%prior = old%prior(permutation)
      if (allocated(old%posterior)) result%posterior = old%posterior(:,permutation)
      if (allocated(old%posterior_unscaled)) result%posterior_unscaled = old%posterior_unscaled(:,permutation)
      if (allocated(old%log_posterior_unscaled)) result%log_posterior_unscaled = old%log_posterior_unscaled(:,permutation)
      if (allocated(old%size)) result%size = old%size(permutation)
      if (allocated(old%beta)) result%beta = old%beta(:,permutation)
      if (allocated(old%multinomial_coef)) result%multinomial_coef = old%multinomial_coef(:,:,permutation)
      if (allocated(old%component_role)) result%component_role = old%component_role(permutation)
      if (allocated(old%sigma)) result%sigma = old%sigma(permutation)
      if (allocated(old%center)) result%center = old%center(:,permutation)
      if (allocated(old%covariance)) result%covariance = old%covariance(:,:,permutation)
      if (allocated(old%factor_loadings)) result%factor_loadings = old%factor_loadings(:,:,permutation)
      if (allocated(old%uniqueness)) result%uniqueness = old%uniqueness(:,permutation)
      if (allocated(old%marginal_variance)) result%marginal_variance = old%marginal_variance(:,permutation)
      if (allocated(old%probability)) result%probability = old%probability(:,permutation)
      if (allocated(old%lambda)) result%lambda = old%lambda(:,permutation)
      if (allocated(old%shape)) result%shape = old%shape(permutation)
      if (allocated(old%rate)) result%rate = old%rate(permutation)
      if (allocated(old%scale)) result%scale = old%scale(permutation)
      if (allocated(old%concomitant_coef)) result%concomitant_coef = old%concomitant_coef(:,permutation)
      if (allocated(old%penalty_lambda)) result%penalty_lambda = old%penalty_lambda(permutation)
      if (allocated(old%smoothing_lambda)) result%smoothing_lambda = old%smoothing_lambda(permutation)
      if (allocated(old%component_effective_df)) result%component_effective_df = old%component_effective_df(permutation)
      if (allocated(old%random_covariance)) result%random_covariance = old%random_covariance(:,:,permutation)
      if (allocated(old%residual_variance)) result%residual_variance = old%residual_variance(permutation)
      if (allocated(old%parameter_design)) result%parameter_design = old%parameter_design(permutation,:)
      if (allocated(old%variance_group)) result%variance_group = old%variance_group(permutation)
      if (allocated(old%cluster)) then
         do i = 1, size(old%cluster)
            result%cluster(i) = inverse(old%cluster(i))
         end do
      end if
      info = 0
   end subroutine flexmix_relabel

   pure subroutine kl_divergence_matrix(object, divergence, eps, overlap)
      real(dp), intent(in) :: object(:,:) !! Nonnegative membership/weight matrix whose columns define discrete distributions.
      real(dp), intent(out) :: divergence(:,:) !! Directional KL-divergence matrix, shape `(k, k)` for `k` input columns.
      real(dp), intent(in), optional :: eps !! Positive floor used before column normalization; default `1e-4` as in FlexMix.
      logical, intent(in), optional :: overlap !! Require at least one jointly above-floor row when true; default true.
      real(dp), allocatable :: p(:,:)
      real(dp) :: floor_value
      logical :: require_overlap, ok
      integer :: i, j
      floor_value = 1.0e-4_dp
      if (present(eps)) floor_value = eps
      require_overlap = .true.
      if (present(overlap)) require_overlap = overlap
      p = max(object, floor_value)
      do j = 1, size(p,2)
         p(:,j) = p(:,j) / sum(p(:,j))
      end do
      divergence = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = 1, size(p,2)
         divergence(i,i) = 0.0_dp
         do j = 1, size(p,2)
            if (i == j) cycle
            ok = any(p(:,i) > floor_value .and. p(:,j) > floor_value)
            if (.not. require_overlap .or. ok) then
               divergence(i,j) = sum(p(:,i) * (log(p(:,i)) - log(p(:,j))))
            end if
         end do
      end do
   end subroutine kl_divergence_matrix

   subroutine kl_divergence_regression(result, x, divergence, info)
      type(flexmix_result), intent(in) :: result !! Fitted Gaussian, Poisson, or binomial regression mixture.
      real(dp), intent(in) :: x(:,:) !! Design matrix on which component conditional distributions are compared.
      real(dp), allocatable, intent(out) :: divergence(:,:) !! Directional component KL-divergence matrix, shape `(k, k)`.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported model kind or prediction failure.
      real(dp), allocatable :: mu(:,:)
      real(dp) :: pk, pl
      integer :: i, j, r
      call flexmix_predict_regression(result, x, mu, info)
      if (info /= 0) then
         allocate(divergence(0,0))
         return
      end if
      allocate(divergence(result%k,result%k))
      divergence = 0.0_dp
      do i = 1, result%k
         do j = 1, result%k
            if (i == j) cycle
            select case (result%model_kind)
            case (flexmix_model_gaussian)
               divergence(i,j) = sum(log(result%sigma(j)) - log(result%sigma(i)) + 0.5_dp * &
                  (-1.0_dp + (result%sigma(i)**2 + (mu(:,i) - mu(:,j))**2) / result%sigma(j)**2))
            case (flexmix_model_poisson)
               divergence(i,j) = sum(mu(:,i) * log(mu(:,i) / mu(:,j)) + mu(:,j) - mu(:,i))
            case (flexmix_model_binomial)
               do r = 1, size(mu,1)
                  pk = clip_probability(mu(r,i))
                  pl = clip_probability(mu(r,j))
                  divergence(i,j) = divergence(i,j) + pk * log(pk / pl) + (1.0_dp - pk) * &
                     log((1.0_dp - pk) / (1.0_dp - pl))
               end do
            case default
               info = 2
               return
            end select
         end do
      end do
      info = 0
   end subroutine kl_divergence_regression

   subroutine kl_divergence_mvnorm(result, divergence, info)
      type(flexmix_result), intent(in) :: result !! Fitted multivariate-normal mixture with centers and covariance matrices.
      real(dp), allocatable, intent(out) :: divergence(:,:) !! Directional Gaussian KL-divergence matrix, shape `(k, k)`.
      integer, intent(out) :: info !! Zero on success; nonzero for wrong model type or non-positive-definite covariance.
      real(dp), allocatable :: inv(:,:), diff(:), work(:,:)
      real(dp) :: logdet_i, logdet_j, trace_value, quad
      integer :: a, b, d, inv_info, r
      if (result%model_kind /= flexmix_model_mvnorm .or. .not. allocated(result%covariance)) then
         allocate(divergence(0,0))
         info = 1
         return
      end if
      d = result%d
      allocate(divergence(result%k,result%k), inv(d,d), diff(d), work(d,d))
      divergence = 0.0_dp
      do a = 1, result%k
         do b = 1, result%k
            if (a == b) cycle
            call inverse_logdet_spd(result%covariance(:,:,b), inv, logdet_j, inv_info)
            if (inv_info /= 0) then
               info = 2
               return
            end if
            call inverse_logdet_spd(result%covariance(:,:,a), work, logdet_i, inv_info)
            if (inv_info /= 0) then
               info = 2
               return
            end if
            work = matmul(inv, result%covariance(:,:,a))
            trace_value = 0.0_dp
            do r = 1, d
               trace_value = trace_value + work(r,r)
            end do
            diff = result%center(:,b) - result%center(:,a)
            quad = dot_product(diff, matmul(inv, diff))
            divergence(a,b) = 0.5_dp * (logdet_j - logdet_i - real(d,dp) + trace_value + quad)
         end do
      end do
      info = 0
   end subroutine kl_divergence_mvnorm

   subroutine step_flexmix_gaussian(x, y, k_values, nrep, step, control)
      real(dp), intent(in) :: x(:,:) !! Gaussian-regression design matrix shared by all candidate fits, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      integer, intent(in) :: k_values(:) !! Candidate component counts to fit, each at least one.
      integer, intent(in) :: nrep !! Number of deterministic initial partitions tried for each candidate component count.
      type(flexmix_step_result), intent(out) :: step !! Candidate log likelihood table and best retained fit for each `k`.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls passed to every candidate fit.
      integer, allocatable :: initial_cluster(:)
      type(flexmix_result) :: fit
      real(dp) :: best
      integer :: i, r, best_r
      step%nrep = max(1,nrep)
      allocate(step%k(size(k_values)), step%logliks(size(k_values),step%nrep), step%models(size(k_values)))
      step%k = k_values
      do i = 1, size(k_values)
         best = -huge(1.0_dp)
         best_r = 1
         do r = 1, step%nrep
            allocate(initial_cluster(size(y)))
            call deterministic_start(y, k_values(i), r, initial_cluster)
            call flexmix_gaussian(x, y, k_values(i), fit, control, initial_cluster=initial_cluster)
            step%logliks(i,r) = fit%loglik
            if (fit%loglik > best .or. r == 1) then
               best = fit%loglik
               best_r = r
               step%models(i) = fit
            end if
            deallocate(initial_cluster)
         end do
         if (best_r < 1) step%models(i) = fit
      end do
   end subroutine step_flexmix_gaussian

   subroutine step_flexmix_unique(step, unique_step)
      type(flexmix_step_result), intent(in) :: step !! Stepwise candidate collection whose duplicate retained component counts.
      type(flexmix_step_result), intent(out) :: unique_step !! New collection retaining the maximum-log-likelihood fit for each.
      integer, allocatable :: retained_k(:), selected(:), order(:)
      integer :: i, j, nmodel, nuniq, pos, tmp

      unique_step%nrep = step%nrep
      if (.not. allocated(step%models) .or. size(step%models) == 0) then
         allocate(unique_step%k(0), unique_step%logliks(0,max(0,step%nrep)), unique_step%models(0))
         return
      end if
      nmodel = size(step%models)
      allocate(retained_k(nmodel), selected(nmodel))
      retained_k = [(step%models(i)%k, i = 1, nmodel)]
      selected = 0
      nuniq = 0
      do i = 1, nmodel
         pos = 0
         do j = 1, nuniq
            if (retained_k(selected(j)) == retained_k(i)) then
               pos = j
               exit
            end if
         end do
         if (pos == 0) then
            nuniq = nuniq + 1
            selected(nuniq) = i
         else if (step%models(i)%loglik > step%models(selected(pos))%loglik) then
            selected(pos) = i
         end if
      end do
      allocate(order(nuniq))
      order = [(i, i = 1, nuniq)]
      do i = 2, nuniq
         tmp = order(i)
         pos = i - 1
         do while (pos >= 1)
            if (retained_k(selected(order(pos))) <= retained_k(selected(tmp))) exit
            order(pos+1) = order(pos)
            pos = pos - 1
         end do
         order(pos+1) = tmp
      end do
      allocate(unique_step%k(nuniq), unique_step%logliks(nuniq,step%nrep), unique_step%models(nuniq))
      do j = 1, nuniq
         i = selected(order(j))
         unique_step%k(j) = step%models(i)%k
         unique_step%models(j) = step%models(i)
         if (allocated(step%logliks) .and. size(step%logliks,1) >= i .and. size(step%logliks,2) == step%nrep) then
            unique_step%logliks(j,:) = step%logliks(i,:)
         else
            unique_step%logliks(j,:) = step%models(i)%loglik
         end if
      end do
   end subroutine step_flexmix_unique

   subroutine get_model(step, criterion, model, info)
      type(flexmix_step_result), intent(in) :: step !! Collection of best fits across candidate component counts.
      character(len=*), intent(in) :: criterion !! Selection criterion: `logLik`, `AIC`, `BIC`, or `ICL`.
      type(flexmix_result), intent(out) :: model !! Candidate minimizing the criterion, or maximizing log likelihood.
      integer, intent(out) :: info !! Zero on success; nonzero for an empty step object or unknown criterion.
      real(dp) :: score, best
      integer :: i, best_i
      if (.not. allocated(step%models) .or. size(step%models) == 0) then
         info = 1
         return
      end if
      best_i = 1
      select case (trim(adjustl(criterion)))
      case ('logLik', 'loglik')
         best = step%models(1)%loglik
         do i = 2, size(step%models)
            if (step%models(i)%loglik > best) then
               best = step%models(i)%loglik
               best_i = i
            end if
         end do
      case ('AIC', 'aic')
         best = flexmix_aic(step%models(1))
         do i = 2, size(step%models)
            score = flexmix_aic(step%models(i))
            if (score < best) then
               best = score
               best_i = i
            end if
         end do
      case ('BIC', 'bic')
         best = flexmix_bic(step%models(1))
         do i = 2, size(step%models)
            score = flexmix_bic(step%models(i))
            if (score < best) then
               best = score
               best_i = i
            end if
         end do
      case ('ICL', 'icl')
         best = flexmix_icl(step%models(1))
         do i = 2, size(step%models)
            score = flexmix_icl(step%models(i))
            if (score < best) then
               best = score
               best_i = i
            end if
         end do
      case default
         info = 2
         return
      end select
      model = step%models(best_i)
      info = 0
   end subroutine get_model

   pure subroutine deterministic_start(signal, k, rep, cluster)
      real(dp), intent(in) :: signal(:) !! Scalar signal whose ranks are partitioned to generate deterministic EM starts.
      integer, intent(in) :: k !! Requested number of starting clusters.
      integer, intent(in) :: rep !! One-based start index controlling alternate deterministic partitions.
      integer, intent(out) :: cluster(:) !! Generated one-based component labels, size equal to `signal`.
      integer, allocatable :: order(:)
      integer :: i, pos, tmp, n, target
      n = size(signal)
      allocate(order(n))
      order = [(i, i = 1, n)]
      do i = 2, n
         tmp = order(i)
         pos = i - 1
         do while (pos >= 1)
            if (signal(order(pos)) <= signal(tmp)) exit
            order(pos+1) = order(pos)
            pos = pos - 1
         end do
         order(pos+1) = tmp
      end do
      if (rep == 1) then
         do i = 1, n
            target = min(k, 1 + (i - 1) * k / max(1,n))
            cluster(order(i)) = target
         end do
      else
         do i = 1, n
            target = 1 + mod((i - 1) * max(1,rep - 1) + rep - 2, k)
            cluster(order(i)) = target
         end do
      end if
   end subroutine deterministic_start

   pure subroutine flexmix_check_result(result, valid)
      type(flexmix_result), intent(in) :: result !! Mixture result whose structural consistency is checked.
      logical, intent(out) :: valid !! True when core dimensions, priors, posteriors, and labels are mutually consistent.
      valid = result%k >= 1 .and. result%n >= 1
      if (.not. valid) return
      if (.not. allocated(result%prior) .or. .not. allocated(result%posterior)) then
         valid = .false.
         return
      end if
      valid = size(result%prior) == result%k .and. size(result%posterior,1) == result%n .and. &
              size(result%posterior,2) == result%k
      if (.not. valid) return
      valid = all(result%prior >= 0.0_dp) .and. abs(sum(result%prior) - 1.0_dp) <= 1.0e-6_dp
      if (allocated(result%cluster)) then
         valid = valid .and. size(result%cluster) == result%n
         if (valid) valid = all(result%cluster >= 1) .and. all(result%cluster <= result%k)
      end if
   end subroutine flexmix_check_result

   subroutine flexmix_remove_component(result, component, info)
      type(flexmix_result), intent(inout) :: result !! Fitted mixture result from which one retained component is removed.
      integer, intent(in) :: component !! One-based retained component index to delete; at least two components must be present.
      integer, intent(out) :: info !! Zero on success; nonzero when deletion is impossible or the index is invalid.
      type(flexmix_result) :: old
      integer, allocatable :: keep(:)
      real(dp) :: rowsum
      integer :: a, i, j, cdf, qdf
      if (result%k <= 1 .or. component < 1 .or. component > result%k) then
         info = 1
         return
      end if
      old = result
      allocate(keep(result%k - 1))
      a = 0
      do j = 1, result%k
         if (j == component) cycle
         a = a + 1
         keep(a) = j
      end do
      result%k = old%k - 1
      result%prior = old%prior(keep)
      result%prior = result%prior / sum(result%prior)
      result%posterior = old%posterior(:,keep)
      do i = 1, result%n
         rowsum = sum(result%posterior(i,:))
         if (rowsum > 0.0_dp) result%posterior(i,:) = result%posterior(i,:) / rowsum
      end do
      if (allocated(old%posterior_unscaled)) result%posterior_unscaled = old%posterior_unscaled(:,keep)
      if (allocated(old%log_posterior_unscaled)) result%log_posterior_unscaled = old%log_posterior_unscaled(:,keep)
      if (allocated(old%beta)) result%beta = old%beta(:,keep)
      if (allocated(old%multinomial_coef)) result%multinomial_coef = old%multinomial_coef(:,:,keep)
      if (allocated(old%component_role)) result%component_role = old%component_role(keep)
      if (allocated(old%sigma)) result%sigma = old%sigma(keep)
      if (allocated(old%center)) result%center = old%center(:,keep)
      if (allocated(old%covariance)) result%covariance = old%covariance(:,:,keep)
      if (allocated(old%factor_loadings)) result%factor_loadings = old%factor_loadings(:,:,keep)
      if (allocated(old%uniqueness)) result%uniqueness = old%uniqueness(:,keep)
      if (allocated(old%marginal_variance)) result%marginal_variance = old%marginal_variance(:,keep)
      if (allocated(old%probability)) result%probability = old%probability(:,keep)
      if (allocated(old%lambda)) result%lambda = old%lambda(:,keep)
      if (allocated(old%shape)) result%shape = old%shape(keep)
      if (allocated(old%rate)) result%rate = old%rate(keep)
      if (allocated(old%scale)) result%scale = old%scale(keep)
      if (allocated(old%concomitant_coef)) result%concomitant_coef = old%concomitant_coef(:,keep)
      if (allocated(old%penalty_lambda)) result%penalty_lambda = old%penalty_lambda(keep)
      if (allocated(old%smoothing_lambda)) result%smoothing_lambda = old%smoothing_lambda(keep)
      if (allocated(old%component_effective_df)) then
         result%component_effective_df = old%component_effective_df(keep)
         result%effective_df = sum(result%component_effective_df)
      end if
      if (allocated(old%random_covariance)) result%random_covariance = old%random_covariance(:,:,keep)
      if (allocated(old%residual_variance)) result%residual_variance = old%residual_variance(keep)
      if (allocated(old%parameter_design)) result%parameter_design = old%parameter_design(keep,:)
      if (allocated(old%variance_group)) result%variance_group = old%variance_group(keep)
      if (allocated(old%cluster)) then
         result%size = old%size(keep)
         result%size = 0
         do i = 1, result%n
            result%cluster(i) = maxloc(result%posterior(i,:), dim=1)
            if (allocated(result%case_weights)) then
               result%size(result%cluster(i)) = result%size(result%cluster(i)) + nint(result%case_weights(i))
            else
               result%size(result%cluster(i)) = result%size(result%cluster(i)) + 1
            end if
         end do
      end if
      cdf = -1
      if (allocated(old%component_role)) then
         if (old%component_role(component) /= flexmix_role_regular) cdf = 0
      end if
      if (cdf < 0) then
         select case (result%model_kind)
         case (flexmix_model_gaussian, flexmix_model_robust_gaussian, flexmix_model_lmc)
            cdf = result%p + 1
         case (flexmix_model_poisson, flexmix_model_binomial, flexmix_model_ziglm_poisson, &
               flexmix_model_ziglm_binomial, flexmix_model_robust_poisson, flexmix_model_conditional_logit)
            cdf = result%p
         case (flexmix_model_gamma_regression)
            cdf = result%p + 1
         case (flexmix_model_multinomial)
            cdf = result%p * max(0,result%nclass - 1)
         case (flexmix_model_lmm, flexmix_model_lmer, flexmix_model_lmmc)
            if (allocated(old%random_covariance)) then
               cdf = result%p + size(old%random_covariance,1) * &
                     (size(old%random_covariance,1) + 1) / 2 + 1
            else
               cdf = result%p + 1
            end if
         case (flexmix_model_mvnorm)
            if (result%diagonal_covariance) then
               cdf = 2 * result%d
            else
               cdf = (3 * result%d + result%d * result%d) / 2
            end if
         case (flexmix_model_factanal)
            cdf = (result%factors + 2) * result%d
         case (flexmix_model_mvbinary, flexmix_model_mvpois)
            cdf = result%d
         case (flexmix_model_mvcombi)
            cdf = result%d + count(.not. result%binary_mask)
         case (flexmix_model_exponential)
            cdf = 1
         case default
            cdf = 2
         end select
      end if
      if (allocated(result%concomitant_coef)) then
         qdf = size(result%concomitant_coef,1)
      else
         qdf = 1
      end if
      result%df = max(0, old%df - cdf - qdf)
      if (allocated(result%component_effective_df)) then
         result%effective_df = sum(result%component_effective_df) + real(qdf*max(0,result%k-1),dp)
         result%df = nint(result%effective_df)
      end if
      if (allocated(result%component_role)) then
         if (all(result%component_role == flexmix_role_regular)) then
            select case (old%model_kind)
            case (flexmix_model_ziglm_poisson, flexmix_model_robust_poisson)
               result%model_kind = flexmix_model_poisson
            case (flexmix_model_ziglm_binomial)
               result%model_kind = flexmix_model_binomial
            case (flexmix_model_robust_gaussian)
               result%model_kind = flexmix_model_gaussian
            end select
         end if
      end if
      info = 0
   end subroutine flexmix_remove_component

   subroutine flexmix_refit_gaussian(old, x, y, result, control, case_weights, group, concomitant_x)
      type(flexmix_result), intent(in) :: old !! Existing Gaussian mixture fit whose posterior memberships initialize the refit.
      real(dp), intent(in) :: x(:,:) !! New or original Gaussian regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! New or original Gaussian response vector, size `n`.
      type(flexmix_result), intent(out) :: result !! Re-estimated Gaussian mixture initialized from `old` when dimensions match.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls for the refit.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: group(:) !! Optional grouping labels, size `n`.
      real(dp), intent(in), optional :: concomitant_x(:,:) !! Optional concomitant prior design matrix, shape `(n, q)`.
      if (allocated(old%posterior)) then
         if (size(old%posterior,1) == size(y)) then
            call flexmix_gaussian(x, y, old%k, result, control, case_weights, group, concomitant_x, &
                                 initial_posterior=old%posterior)
         else
            call flexmix_gaussian(x, y, max(1,old%k), result, control, case_weights, group, concomitant_x)
         end if
      else
         call flexmix_gaussian(x, y, max(1,old%k), result, control, case_weights, group, concomitant_x)
      end if
   end subroutine flexmix_refit_gaussian


   subroutine flexmix_make_distribution(family, beta, prior, result, info, sigma, shape)
      character(len=*), intent(in) :: family !! Regression family name: `gaussian`, `poisson`, `binomial`, or `Gamma`.
      real(dp), intent(in) :: beta(:,:) !! Component regression coefficients, shape `(p,k)`.
      real(dp), intent(in) :: prior(:) !! Nonnegative component prior masses, size `k`; normalized internally.
      type(flexmix_result), intent(out) :: result !! Fixed numerical mixture for prediction and parameter access.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, family, or required parameters.
      real(dp), intent(in), optional :: sigma(:) !! Gaussian component standard deviations, size `k`, all strictly positive.
      real(dp), intent(in), optional :: shape(:) !! Gamma component shape parameters, size `k`, all strictly positive.
      integer :: k

      k = size(beta,2)
      result%status = 0
      if (k < 1 .or. size(prior) /= k .or. any(prior < 0.0_dp) .or. sum(prior) <= 0.0_dp) then
         result%status = -1
         info = 1
         return
      end if
      result%n = 0
      result%p = size(beta,1)
      result%d = 1
      result%k = k
      result%k0 = k
      result%prior = prior / sum(prior)
      result%beta = beta
      allocate(result%component_role(k))
      result%component_role = flexmix_role_regular
      select case (trim(adjustl(family)))
      case ('gaussian')
         if (.not. present(sigma)) then
            result%status = -2
            info = 2
            return
         end if
         if (size(sigma) /= k .or. any(sigma <= 0.0_dp)) then
            result%status = -2
            info = 2
            return
         end if
         result%model_kind = flexmix_model_gaussian
         result%sigma = sigma
         result%df = k * (result%p + 1) + k - 1
      case ('poisson')
         result%model_kind = flexmix_model_poisson
         result%df = k * result%p + k - 1
      case ('binomial')
         result%model_kind = flexmix_model_binomial
         result%df = k * result%p + k - 1
      case ('Gamma', 'gamma')
         if (.not. present(shape)) then
            result%status = -3
            info = 3
            return
         end if
         if (size(shape) /= k .or. any(shape <= 0.0_dp)) then
            result%status = -3
            info = 3
            return
         end if
         result%model_kind = flexmix_model_gamma_regression
         result%shape = shape
         result%df = k * (result%p + 1) + k - 1
      case default
         result%status = -4
         result%model_kind = 0
         info = 4
         return
      end select
      info = 0
   end subroutine flexmix_make_distribution

   subroutine flexmix_simulate_regression(result, x, y, cluster, info)
      type(flexmix_result), intent(in) :: result !! Fitted regression mixture supplying component priors and coefficients.
      real(dp), intent(in) :: x(:,:) !! Predictor design matrix for simulated responses, shape `(n, p)`.
      real(dp), intent(out) :: y(:) !! Simulated response vector, size `n`; binomial output is Bernoulli zero/one.
      integer, intent(out) :: cluster(:) !! Simulated one-based component label for each row, size `n`.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported model type or incompatible dimensions.
      real(dp) :: u, mu, p, z
      integer :: i, j
      if (.not. allocated(result%beta) .or. size(x,2) /= result%p .or. size(y) /= size(x,1) .or. &
          size(cluster) /= size(x,1)) then
         info = 1
         return
      end if
      do i = 1, size(x,1)
         call random_number(u)
         cluster(i) = result%k
         p = 0.0_dp
         do j = 1, result%k
            p = p + result%prior(j)
            if (u <= p) then
               cluster(i) = j
               exit
            end if
         end do
         mu = dot_product(x(i,:), result%beta(:,cluster(i)))
         select case (result%model_kind)
         case (flexmix_model_gaussian)
            call standard_normal_random(z)
            y(i) = mu + result%sigma(cluster(i)) * z
         case (flexmix_model_poisson)
            call poisson_random(exp(max(-30.0_dp,min(30.0_dp,mu))), y(i))
         case (flexmix_model_binomial)
            p = 1.0_dp / (1.0_dp + exp(-max(-30.0_dp,min(30.0_dp,mu))))
            call random_number(u)
            if (u < p) then
               y(i) = 1.0_dp
            else
               y(i) = 0.0_dp
            end if
         case default
            info = 2
            return
         end select
      end do
      info = 0
   end subroutine flexmix_simulate_regression

   subroutine ex_linear(beta, n_per_component, family, x, y, cluster, info, sigma)
      real(dp), intent(in) :: beta(:,:) !! Regression coefficient matrix, shape `(d+1, k)`, with intercept in row one.
      integer, intent(in) :: n_per_component(:) !! Number of observations generated from each component, size `k`.
      character(len=*), intent(in) :: family !! Response family, `gaussian` or `poisson`.
      real(dp), allocatable, intent(out) :: x(:,:) !! Generated uniform-on-zero-to-one predictors, shape `(sum(n), d)`.
      real(dp), allocatable, intent(out) :: y(:) !! Generated response vector, size `sum(n)`.
      integer, allocatable, intent(out) :: cluster(:) !! Known one-based generating component labels, size `sum(n)`.
      integer, intent(out) :: info !! Zero on success; nonzero for inconsistent dimensions or unsupported family.
      real(dp), intent(in), optional :: sigma(:) !! Optional Gaussian component standard deviations; defaults to one.
      real(dp), allocatable :: sd(:)
      real(dp) :: eta, z
      integer :: d, i, j, k, n, pos
      k = size(beta,2)
      d = size(beta,1) - 1
      if (size(n_per_component) /= k .or. any(n_per_component < 0) .or. d < 0) then
         allocate(x(0,0), y(0), cluster(0))
         info = 1
         return
      end if
      allocate(sd(k))
      sd = 1.0_dp
      if (present(sigma)) then
         if (size(sigma) /= k .or. any(sigma <= 0.0_dp)) then
            allocate(x(0,0), y(0), cluster(0))
            info = 2
            return
         end if
         sd = sigma
      end if
      n = sum(n_per_component)
      allocate(x(n,d), y(n), cluster(n))
      if (d > 0) call random_number(x)
      pos = 0
      do j = 1, k
         do i = 1, n_per_component(j)
            pos = pos + 1
            cluster(pos) = j
            eta = beta(1,j)
            if (d > 0) eta = eta + dot_product(x(pos,:), beta(2:,j))
            select case (trim(adjustl(family)))
            case ('gaussian')
               call standard_normal_random(z)
               y(pos) = eta + sd(j) * z
            case ('poisson')
               call poisson_random(exp(max(-30.0_dp,min(30.0_dp,eta))), y(pos))
            case default
               y = 0.0_dp
               info = 3
               return
            end select
         end do
      end do
      info = 0
   end subroutine ex_linear


   subroutine ex_npreg(n, x, yn, yp, yb, class_label, id1, id2, info)
      integer, intent(in) :: n !! Per-class sample size; must be positive and even as in the upstream example generator.
      real(dp), allocatable, intent(out) :: x(:) !! Uniform predictors on `[0,10]`, size `2*n`.
      real(dp), allocatable, intent(out) :: yn(:) !! Nonlinear Gaussian responses, size `2*n`.
      real(dp), allocatable, intent(out) :: yp(:) !! Poisson responses, size `2*n`, represented in the package real kind.
      real(dp), allocatable, intent(out) :: yb(:) !! Bernoulli responses, size `2*n`, represented as zero/one reals.
      integer, allocatable, intent(out) :: class_label(:) !! Known generating class labels, first `n` equal to one and last `n`
      integer, allocatable, intent(out) :: id1(:) !! Upstream `id1` grouping labels, pairing consecutive rows across the `2*n`
      integer, allocatable, intent(out) :: id2(:) !! Upstream `id2` grouping labels, grouping consecutive blocks of four rows.
      integer, intent(out) :: info !! Zero on success; nonzero when `n` is not a positive even integer.
      real(dp), allocatable :: noise(:)
      real(dp) :: mean_value, probability, u
      integer :: i

      if (n <= 0 .or. mod(n,2) /= 0) then
         allocate(x(0), yn(0), yp(0), yb(0), class_label(0), id1(0), id2(0))
         info = 1
         return
      end if
      allocate(x(2*n), yn(2*n), yp(2*n), yb(2*n), class_label(2*n), id1(2*n), id2(2*n), noise(n))
      call random_number(x)
      x = 10.0_dp * x
      do i = 1, n
         call standard_normal_random(noise(i))
      end do
      do i = 1, n
         yn(i) = 5.0_dp * x(i) + 3.0_dp * noise(i)
         yn(n+i) = 40.0_dp - (x(n+i) - 5.0_dp)**2 + 3.0_dp * noise(i)
         mean_value = exp(2.0_dp - 0.2_dp * x(i))
         call poisson_random(mean_value, yp(i))
         mean_value = exp(1.0_dp + 0.1_dp * x(n+i))
         call poisson_random(mean_value, yp(n+i))
         probability = 1.0_dp / (1.0_dp + exp(-(x(i) - 5.0_dp)))
         call random_number(u)
         yb(i) = merge(1.0_dp, 0.0_dp, u < probability)
         probability = 1.0_dp / (1.0_dp + exp(-(5.0_dp - x(n+i))))
         call random_number(u)
         yb(n+i) = merge(1.0_dp, 0.0_dp, u < probability)
      end do
      class_label(1:n) = 1
      class_label(n+1:2*n) = 2
      do i = 1, 2*n
         id1(i) = (i + 1) / 2
         id2(i) = (i + 3) / 4
      end do
      info = 0
   end subroutine ex_npreg

   subroutine ex_nclus(n, y, cluster, info)
      integer, intent(in) :: n !! Base cluster size; must be positive and even so the third cluster has exactly `1.5*n` rows.
      real(dp), allocatable, intent(out) :: y(:,:) !! Two-dimensional Gaussian-cluster sample, shape `(11*n/2,2)`.
      integer, allocatable, intent(out) :: cluster(:) !! Known generating cluster labels one through four for each returned row.
      integer, intent(out) :: info !! Zero on success; nonzero when `n` is not a positive even integer.
      real(dp) :: z1, z2
      integer :: i, n3, total, pos

      if (n <= 0 .or. mod(n,2) /= 0) then
         allocate(y(0,0), cluster(0))
         info = 1
         return
      end if
      n3 = 3*n/2
      total = n + n + n3 + 2*n
      allocate(y(total,2), cluster(total))
      pos = 0
      do i = 1, n
         call standard_normal_random(z1)
         call standard_normal_random(z2)
         pos = pos + 1
         y(pos,:) = [z1, z2]
         cluster(pos) = 1
      end do
      do i = 1, n
         call standard_normal_random(z1)
         call standard_normal_random(z2)
         pos = pos + 1
         y(pos,:) = [8.0_dp + z1, sqrt(2.0_dp) * z2]
         cluster(pos) = 2
      end do
      do i = 1, n3
         call standard_normal_random(z1)
         call standard_normal_random(z2)
         pos = pos + 1
         y(pos,:) = [-2.0_dp + sqrt(2.0_dp) * z1, 6.0_dp + z2]
         cluster(pos) = 3
      end do
      do i = 1, 2*n
         call standard_normal_random(z1)
         call standard_normal_random(z2)
         pos = pos + 1
         y(pos,1) = 4.0_dp + z1
         y(pos,2) = 4.0_dp + 0.9_dp * z1 + sqrt(0.19_dp) * z2
         cluster(pos) = 4
      end do
      info = 0
   end subroutine ex_nclus

   subroutine standard_normal_random(z)
      real(dp), intent(out) :: z !! One standard-normal random variate generated by the Box-Muller transform.
      real(dp) :: u1, u2
      real(dp), parameter :: pi_local = acos(-1.0_dp)
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi_local * u2)
   end subroutine standard_normal_random

   subroutine poisson_random(lambda_value, value)
      real(dp), intent(in) :: lambda_value !! Nonnegative Poisson mean for one simulated count.
      real(dp), intent(out) :: value !! Simulated nonnegative integer count represented in the package real kind.
      real(dp) :: limit, product_value, u
      integer :: count_value
      if (lambda_value <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (lambda_value > 30.0_dp) then
         call standard_normal_random(u)
         value = real(max(0,nint(lambda_value + sqrt(lambda_value) * u)),dp)
         return
      end if
      limit = exp(-lambda_value)
      product_value = 1.0_dp
      count_value = -1
      do while (product_value > limit)
         count_value = count_value + 1
         call random_number(u)
         product_value = product_value * max(u, tiny(1.0_dp))
      end do
      value = real(count_value,dp)
   end subroutine poisson_random

end module flexmix_api
