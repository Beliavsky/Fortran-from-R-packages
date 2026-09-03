! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Public Fortran facade for computational algorithms translated from mice 3.19.0.
module mice
    use r_kinds, only : dp
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_singular, &
                            mice_no_observed, mice_not_converged
    use mice_rng, only : mice_rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq
    use mice_native, only : legendre_basis
    use mice_matching, only : matchindex, matcher, pmm_match_value
    use mice_regression, only : mice_norm_draw, estimice, norm_draw, impute_norm, impute_norm_nob, &
                                impute_norm_predict, impute_norm_boot
    use mice_impute_continuous, only : impute_mean, impute_sample, impute_pmm, impute_random_indicator, &
                                       impute_quadratic
    use mice_midastouch, only : impute_midastouch
    use mice_mpmm, only : impute_mpmm
    use mice_categorical, only : logistic_fit, impute_logreg, impute_logreg_boot, multinomial_fit, &
                                 impute_polyreg
    use mice_polr, only : proportional_odds_fit, proportional_odds_probabilities, impute_polr
    use mice_lda, only : impute_lda
    use mice_twolevel, only : impute_2lonly_mean, impute_2lonly_norm, impute_2lonly_pmm
    use mice_2l_norm, only : mice_2l_norm_state, impute_2l_norm
    use mice_mnar, only : impute_mnar_norm, impute_mnar_logreg
    use mice_fcs, only : mice_fcs_result, mice_fcs_impute, mice_method_skip, mice_method_mean, &
                         mice_method_sample, mice_method_norm, mice_method_norm_nob, mice_method_norm_predict, &
                         mice_method_norm_boot, mice_method_pmm, mice_method_logreg, mice_method_logreg_boot, &
                         mice_method_polyreg, mice_method_polr, mice_method_midastouch, mice_method_lda
    use mice_pooling, only : pool_scalar_result, pool_vector_result, d3_result, barnard_rubin, pool_scalar, &
                             pool_vector, pooled_wald, d3_from_deviances
    use mice_diagnostics, only : md_pairs_result, flux_result, md_pattern_result, missing_mask, md_pairs, &
                                 md_pairs_from_mask, flux, md_pattern, quickpred
    use mice_ampute, only : ampute_right, ampute_left, ampute_mid, ampute_tail, ampute_mcar, &
                            ampute_continuous, ampute_discrete
    implicit none
    private

    public :: dp
    public :: mice_ok, mice_invalid_argument, mice_invalid_shape, mice_singular, mice_no_observed, mice_not_converged
    public :: mice_rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq
    public :: legendre_basis, matchindex, matcher, pmm_match_value
    public :: mice_norm_draw, estimice, norm_draw, impute_norm, impute_norm_nob, impute_norm_predict, impute_norm_boot
    public :: impute_mean, impute_sample, impute_pmm, impute_random_indicator, impute_quadratic
    public :: impute_midastouch, impute_mpmm
    public :: logistic_fit, impute_logreg, impute_logreg_boot, multinomial_fit, impute_polyreg, impute_polr
    public :: proportional_odds_fit, proportional_odds_probabilities, impute_lda
    public :: impute_2lonly_mean, impute_2lonly_norm, impute_2lonly_pmm
    public :: mice_2l_norm_state, impute_2l_norm, impute_mnar_norm, impute_mnar_logreg
    public :: mice_fcs_result, mice_fcs_impute
    public :: mice_method_skip, mice_method_mean, mice_method_sample, mice_method_norm, mice_method_norm_nob
    public :: mice_method_norm_predict, mice_method_norm_boot, mice_method_pmm, mice_method_logreg
    public :: mice_method_logreg_boot, mice_method_polyreg, mice_method_polr, mice_method_midastouch, mice_method_lda
    public :: pool_scalar_result, pool_vector_result, d3_result, barnard_rubin, pool_scalar, pool_vector, pooled_wald
    public :: d3_from_deviances
    public :: md_pairs_result, flux_result, md_pattern_result, missing_mask, md_pairs, md_pairs_from_mask, flux
    public :: md_pattern, quickpred
    public :: ampute_right, ampute_left, ampute_mid, ampute_tail, ampute_mcar, ampute_continuous, ampute_discrete

end module mice
