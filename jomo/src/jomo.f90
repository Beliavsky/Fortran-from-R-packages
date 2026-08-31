! Public facade for the modern Fortran translation of the jomo computational core.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq, rng_student_t
   use jomo_single_level, only : jomo1_result, jomo1_mixed_mcmc, jomo1con_mcmc, jomo1cat_mcmc
   use jomo_multilevel, only : jomo1ran_result, jomo1ran_mixed_mcmc, jomo1rancon_mcmc
   use jomo_heteroscedastic, only : jomo1ranhr_result, jomo1ranhr_mixed_mcmc, jomo1ranconhr_mcmc
   use jomo_twolevel, only : jomo2_result, jomo2_mixed_mcmc, jomo2con_mcmc
   use jomo_twolevel_heteroscedastic, only : jomo2hr_result, jomo2hr_mixed_mcmc, jomo2conhr_mcmc
   use jomo_substantive, only : normal_cdf, linear_loglik, binary_probit_loglik, ordinal_probit_loglik
   use jomo_substantive, only : cox_partial_loglik_ordered, sample_binary_probit_latent
   use jomo_substantive, only : sample_ordinal_probit_latent, update_ordinal_thresholds
   use jomo_substantive, only : sample_gaussian_coefficients, sample_linear_variance, cox_coordinate_newton
   use jomo_smc, only : smc_linear, smc_binary_probit, smc_cox, smc_ordinal_probit
   use jomo_smc, only : smc_l1_continuous, smc_l1_categorical, smc_l2_continuous, smc_l2_categorical
   use jomo_smc, only : smc_factor_spec, smc_term_spec, smc_design_spec, smc_substantive_model, smc_sweep_stats
   use jomo_smc, only : smc_design_columns, smc_build_design, smc_substantive_loglik
   use jomo_smc, only : smc_level1_sweep, smc_level2_sweep, smc_update_substantive, smc_compatible_iteration
   implicit none
   private

   public :: dp, i8
   public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq, rng_student_t
   public :: jomo1_result, jomo1_mixed_mcmc, jomo1con_mcmc, jomo1cat_mcmc
   public :: jomo1ran_result, jomo1ran_mixed_mcmc, jomo1rancon_mcmc
   public :: jomo1ranhr_result, jomo1ranhr_mixed_mcmc, jomo1ranconhr_mcmc
   public :: jomo2_result, jomo2_mixed_mcmc, jomo2con_mcmc
   public :: jomo2hr_result, jomo2hr_mixed_mcmc, jomo2conhr_mcmc
   public :: normal_cdf, linear_loglik, binary_probit_loglik, ordinal_probit_loglik
   public :: cox_partial_loglik_ordered, sample_binary_probit_latent
   public :: sample_ordinal_probit_latent, update_ordinal_thresholds
   public :: sample_gaussian_coefficients, sample_linear_variance, cox_coordinate_newton
   public :: smc_linear, smc_binary_probit, smc_cox, smc_ordinal_probit
   public :: smc_l1_continuous, smc_l1_categorical, smc_l2_continuous, smc_l2_categorical
   public :: smc_factor_spec, smc_term_spec, smc_design_spec, smc_substantive_model, smc_sweep_stats
   public :: smc_design_columns, smc_build_design, smc_substantive_loglik
   public :: smc_level1_sweep, smc_level2_sweep, smc_update_substantive, smc_compatible_iteration
end module jomo
