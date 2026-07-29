! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff
   use fracdiff_kinds, only : dp, pi_dp
   use fracdiff_status
   use fracdiff_types
   use fracdiff_difference, only : diffseries, diffseries_direct, fractional_weights
   use fracdiff_estimators, only : fd_gph, fd_sperio
   use fracdiff_simulation_mod, only : fracdiff_sim, fractional_arma_filter_simulation
   use fracdiff_model_api, only : fracdiff_fit, fracdiff_var, fracdiff_coefficients, &
                                  fracdiff_confint, summarize_fracdiff, fracdiff_aic, fracdiff_bic
   use fracdiff_filter, only : haslett_raftery_filter, arma_residuals, &
                               arma_residual_jacobian, conditional_arma_series_residuals
   use fracdiff_polynomial, only : polynomial_multiply, polynomial_roots, &
                                   minimum_ar_root_modulus
   implicit none
   public

end module fracdiff
