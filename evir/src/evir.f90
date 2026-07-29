! SPDX-License-Identifier: GPL-2.0-or-later
module evir
    use evir_kinds, only : dp
    use evir_types
    use evir_math, only : seed_rng, random_uniform, random_normal
    use evir_distributions
    use evir_data, only : findthresh => find_threshold, block_maxima, block_maxima_groups, decluster
    use evir_fitting, only : gev => fit_gev, gumbel => fit_gumbel, gpd => fit_gpd, pot => fit_pot, &
        rlevel_gev => gev_return_level, rlevel_gev_profile => gev_return_level_profile, &
        gpd_q => gpd_quantile_estimate, gpd_q_wald => gpd_quantile_wald, &
        gpd_q_profile => gpd_quantile_profile, gpd_sfall => gpd_shortfall_estimate, &
        gpd_sfall_profile => gpd_shortfall_profile, riskmeasures => risk_measures, &
        gev_negloglik_value, gpd_negloglik_value, pot_negloglik_value
    use evir_eda, only : emplot => empirical_tail, meplot => mean_excess, &
        qplot => quantile_plot_data, records => records_analysis, hill => hill_estimates, &
        exindex => extremal_index, shape => shape_stability, quant => quantile_stability, &
        tailplot => tail_curve
    use evir_bivariate, only : gpdbiv => fit_gpdbiv, interpret_gpdbiv, &
        logistic_exponent, bivariate_cdf, bivariate_survivor
    implicit none
    public
end module evir
