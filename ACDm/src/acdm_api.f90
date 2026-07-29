! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_api
  use acdm_fit, only : acdFit => acd_fit_model, acd_fit_options, acd_fit_result, &
                       forecast_acd, acd_bootstrap_se, acd_score_matrix
  use acdm_models, only : sim_ACD => simulate_acd, acd_order, model_code, &
                          model_name, filter_acd, acd_loglik
  use acdm_data, only : computeDurations => compute_durations, &
                        diurnalAdj => diurnal_adjust, duration_result, &
                        diurnal_result, DURATION_TRADE, DURATION_PRICE, &
                        DURATION_VOLUME, DIURNAL_CUBIC_SPLINE, &
                        DIURNAL_SMOOTH_SPLINE, DIURNAL_SUPER_SMOOTHER, &
                        DIURNAL_FFF
  use acdm_diagnostics, only : standardizeResi => standardize_residuals, &
                               acf_acd, resiDensityAcd => residual_density_acd, &
                               qqplotAcd => qqplot_acd, &
                               testRmACD => test_rm_acd, &
                               testSTACD => test_st_acd, &
                               testTVACD => test_tv_acd, &
                               lm_test_result, acf_result, density_result, &
                               qq_result, summary_result
  use acdm_distributions, only : dburr, pburr, qburr, rburr, &
       burrExpectation => burr_expectation, dgenf, pgenf, qgenf, rgenf, &
       genfHazard => genf_hazard, dgengamma, pgengamma, qgengamma, &
       rgengamma, gengammaHazard => gengamma_hazard, dqweibull, &
       pqweibull, qqweibull, rqweibull, &
       qweibullExpectation => qweibull_expectation, &
       qweibullHazard => qweibull_hazard, dmixqwe, pmixqwe, qmixqwe, &
       rmixqwe, mixqweHazard => mixqwe_hazard, dmixqww, pmixqww, &
       qmixqww, rmixqww, mixqwwHazard => mixqww_hazard, dmixinvgauss, &
       pmixinvgauss, qmixinvgauss, rmixinvgauss, &
       mixinvgaussHazard => mixinvgauss_hazard
  use acdm_profiles, only : plotHazard => hazard_diagnostics, &
                            plotLL => likelihood_profile, hazard_result, &
                            likelihood_profile_result
  implicit none
  public
end module acdm_api
