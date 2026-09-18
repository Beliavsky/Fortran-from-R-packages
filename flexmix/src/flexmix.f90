! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix
   use flexmix_kinds, only : dp
   use flexmix_types
   use flexmix_api
   use flexmix_components, only : fit_gaussian_regression_component, fit_poisson_regression_component
   use flexmix_components, only : fit_binomial_regression_component, fit_mvnormal_component
   use flexmix_components, only : fit_conditional_logit_component, conditional_logit_log_density
   use flexmix_components, only : fit_factor_analysis_component
   use flexmix_components, only : fit_mvbinary_component, fit_mvpois_component, fit_mvcombi_component
   use flexmix_components, only : fit_lognormal_component, fit_exponential_component
   use flexmix_components, only : fit_inverse_gaussian_component, fit_gamma_component, fit_weibull_component
   use flexmix_components, only : gaussian_regression_log_density, poisson_regression_log_density
   use flexmix_components, only : binomial_regression_log_density, mvnormal_log_density
   use flexmix_components, only : mvbinary_log_density, mvpois_log_density, mvcombi_log_density
   use flexmix_components, only : lognormal_log_density, exponential_log_density
   use flexmix_components, only : inverse_gaussian_log_density, gamma_log_density, weibull_log_density
   use flexmix_em_utils, only : initialize_posteriors, classify_posteriors, group_log_densities, group_first_mask
   use flexmix_em_utils, only : fit_constant_prior, fit_multinomial_prior, multinomial_prior
   use flexmix_refit
   use flexmix_penalized
   use flexmix_smooth
   use flexmix_mixed
   use flexmix_fixed
   use flexmix_bootstrap
   implicit none
   public
end module flexmix
