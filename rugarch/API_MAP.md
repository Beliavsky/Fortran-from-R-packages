# API coverage map

This table maps important R `rugarch` numerical entry points to this
experimental Fortran project. The translation exposes typed numerical APIs
rather than reproducing R classes and dispatch.

| R API | Fortran API | Status |
|---|---|---|
| `ugarchspec` | `garch_spec`, `make_garch_spec` | Typed Fortran specification |
| `ugarchfilter` | `garch_filter`, `realgarch_filter`, `filter_garch_extended` | Main and advanced recursions; external-regressor wrapper included |
| `ugarchfit` | `fit_model`, specialized `fit_*`, `fit_garch_extended` | Direct fitting plus regressors, ARCH-in-mean, targeting, and covariance output |
| `ugarchsim`, `ugarchpath` | `simulate_garch`, `simulate_realgarch` | Implemented for supported recursions |
| `ugarchforecast` | `forecast_volatility`, `garch_bootstrap_forecast` | Analytic volatility and conditional bootstrap forecast densities |
| `ugarchboot` | `garch_bootstrap_forecast` | Raw, kernel, and semi-parametric sampling; partial and full modes |
| `ugarchdistribution` | `garch_parametric_distribution` | Parametric simulate/refit distribution |
| `uncvariance`, `uncmean`, `halflife` | `unconditional_variance`, `unconditional_mean`, `volatility_half_life` | Implemented/approximate by model |
| `persistence` | `true_persistence` | Implemented |
| `newsimpact` | `news_impact` | Implemented for common models |
| mean/variance external regressors | `fit_garch_extended`, `filter_garch_extended` | Experimental dependency-free numerical equivalent |
| variance targeting | `fit_garch_extended(...,variance_targeting=.true.)` | Implemented for supported model families |
| fit covariance/robust SE | `attach_garch_covariance` and automatic extended-fit attachment | Numerical Hessian and Newey-West sandwich covariance |
| `ddist`, `pdist`, `qdist`, `rdist` | `distribution_pdf`, `distribution_cdf`, `distribution_quantile`, `random_innovation` | norm/snorm/std/sstd/ged/sged/jsu/NIG/GHYP/GHST |
| `fitdist` | `fit_distribution` | Implemented for supported distributions |
| `dskewness`, `dkurtosis` | `distribution_skewness`, `distribution_excess_kurtosis` | Analytic where available, otherwise numerical quadrature |
| `ghyptransform` | `ghyp_transform` | Standardized-to-raw GH parameters |
| GHYP functions | `dsgh`, `psgh`, `qsgh`, `rsgh` | Native standardized GH implementation |
| NIG functions | `dsnig`, `psnig`, `qsnig`, `rsnig` | Native standardized NIG implementation |
| GH skew-t functions | `dsghst`, `psghst`, `qsghst`, `rsghst` | Native standardized GH skew-t implementation |
| FIGARCH model | `figarch_weights`, `garch_filter`, `fit_figarch` | Arbitrary p/q, truncated infinite-ARCH recursion |
| component-GARCH model | `garch_filter`, `fit_csgarch` | Arbitrary p/q direct likelihood fit |
| realGARCH model | `realgarch_filter`, `realgarch_log_likelihood`, `simulate_realgarch`, `fit_realgarch` | Joint return/measurement fit |
| fGARCH model | `configure_fgarch_submodel`, `garch_filter`, `fit_fgarch` | Eight rugarch Hentschel submodels |
| `arfimaspec` | `arfima_spec`, `make_arfima_spec` | Implemented |
| `arfimafilter` | `arma_residuals`, `fractional_difference` | Numeric filtering implemented |
| `arfimasim`, `arfimapath` | `simulate_arfima` | Implemented |
| `arfimafit`, `autoarfima` | `fit_arfima`, `gph_estimate_d`, `auto_arfima` | Approximate AIC/BIC order search |
| `arfimaforecast` | `forecast_arfima`, `arfima_bootstrap_forecast` | Conditional and residual-bootstrap forecasts |
| `arfimaroll` | `rolling_arfima_forecast` | Serial rolling implementation |
| `arfimadistribution` | `arfima_parametric_distribution` | Parametric simulate/refit distribution |
| `arfimacv` | `arfima_cross_validation` | Expanding-window p/q comparison |
| ARFIMA multi methods | `multifit_arfima`, `multiforecast_arfima` | Sequential numeric wrappers |
| GARCH multi methods | `multifit_garch`, `multiforecast_garch` | Sequential numeric wrappers |
| `ugarchroll` | `rolling_garch_forecast` | Serial one-step rolling forecasts |
| `VaRTest` | `var_test` | Kupiec UC and Christoffersen CC |
| `ESTest` | `es_test_p_value` | Implemented nonparametric version |
| `DACTest` | `directional_accuracy_test` | Pesaran-Timmermann version |
| `BerkowitzTest` | `berkowitz_test` | AR(1) likelihood-ratio version |
| `VaRloss` | `var_loss`, `quantile_loss` | Implemented |
| `mcsTest` | `mcs_test` | Range and semi-quadratic MCS elimination |
| `GMMTest`, `HLTest`, `VaRDurTest` | `gmm_test`, `hong_li_test`, `var_duration_test` | Implemented |
| `nyblom`, `signbias`, `gof` | `nyblom_test`, `sign_bias_test`, `adjusted_pearson_gof` | Implemented from numeric inputs |
| weighted residual/ARCH tests | `weighted_box_test`, `arch_lm_test` | Implemented |
| forecast performance | `forecast_performance` and loss functions | Implemented |
| `move`, `generatefwd`, `ftseq` numeric core | numeric index utilities | Numeric equivalents without R date classes |
| R S4/S3 methods and plotting | — | Intentionally omitted |
| parallel rolling | — | Intentionally omitted |
| resumable/checkpoint rolling | — | Intentionally omitted |
