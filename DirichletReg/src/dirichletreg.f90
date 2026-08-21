! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg
  use dirichletreg_kinds, only : dp
  use dirichletreg_types, only : design_block, dirichletreg_model
  use dirichletreg_rng, only : seed_rng
  use dirichletreg_distribution, only : ddirichlet, ddirichlet_log, rdirichlet
  use dirichletreg_data, only : prepare_composition, prepare_beta, collapse_subcomposition
  use dirichletreg_common, only : common_loglik_score, common_loglik_score_hessian, common_predict
  use dirichletreg_alternative, only : alternative_loglik_score, alternative_loglik_score_hessian, alternative_predict
  use dirichletreg_fit, only : fit_common, fit_alternative
  use dirichletreg_inference, only : standardized_residuals, raw_residuals, composite_residuals, &
       wald_confint, coefficient_tests, likelihood_ratio_test
  use dirichletreg_geometry, only : to_ternary, to_quaternary
  implicit none
  public
end module dirichletreg
