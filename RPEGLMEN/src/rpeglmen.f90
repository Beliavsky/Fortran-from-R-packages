! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

module rpeglmen
  use rpeglmen_kinds, only : dp
  use rpeglmen_types, only : enet_options, fit_result, path_result, &
    rpe_success, rpe_invalid_input, rpe_no_convergence, rpe_line_search_failure, &
    rpe_numerical_failure, model_exponential, model_gamma
  use rpeglmen_penalty, only : prox_l1, prox_en, regularizer_en
  use rpeglmen_likelihood, only : exp_negative_log_likelihood, &
    grad_exp_negative_log_likelihood, gamma_negative_log_likelihood, &
    grad_gamma_negative_log_likelihood, predict_mean, mse_value, mse_gradient
  use rpeglmen_solver, only : fit_fixed_model, glmnet_exp_fixed, glm_gamma_net_fixed
  use rpeglmen_gamma_mle, only : fit_glm_gamma_mle
  use rpeglmen_cv, only : compute_lambda_max, generate_lambda_grid, &
    fit_regularization_path, glmnet_exp, fit_glm_gamma_net
  use rpeglmen_compat, only : fit_glmGammaNet, cv_glmGammaNet, glmGammaNet, &
    fitGlmCv, ExpNegativeLogLikelihood_cpp, GradExpNegativeLogLikelihood_cpp, &
    ProxGradDescent_cpp
  implicit none
  private

  public :: dp
  public :: enet_options, fit_result, path_result
  public :: rpe_success, rpe_invalid_input, rpe_no_convergence
  public :: rpe_line_search_failure, rpe_numerical_failure
  public :: model_exponential, model_gamma
  public :: prox_l1, prox_en, regularizer_en
  public :: exp_negative_log_likelihood, grad_exp_negative_log_likelihood
  public :: gamma_negative_log_likelihood, grad_gamma_negative_log_likelihood
  public :: predict_mean, mse_value, mse_gradient
  public :: fit_fixed_model, glmnet_exp_fixed, glm_gamma_net_fixed
  public :: fit_glm_gamma_mle, compute_lambda_max, generate_lambda_grid
  public :: fit_regularization_path, glmnet_exp, fit_glm_gamma_net
  public :: fit_glmGammaNet, cv_glmGammaNet, glmGammaNet, fitGlmCv
  public :: ExpNegativeLogLikelihood_cpp, GradExpNegativeLogLikelihood_cpp
  public :: ProxGradDescent_cpp

end module rpeglmen
