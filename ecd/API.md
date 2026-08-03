# API overview

All public numerical procedures use double precision from `ecd_kinds`.
Status values are:

- `ecd_ok = 0`
- `ecd_invalid = 1`
- `ecd_no_convergence = 2`

## `ecd_api`

Umbrella module that re-exports the complete public API.

## `ecd_core`

Types:

- `ecd_model`
- `ecd_stats_type`
- `ellipticity_result`

Construction and parameter transformations:

- `ecd_new`, `ecd_cusp_new`, `ecd_polar_new`, `ecd_initialize`
- `ecd_cusp_a2r`, `ecd_cusp_r2a`
- `ecd_adj_gamma`, `ecd_adj2gamma`
- `ecd_discriminant`, `ecd_j_invariant`, `ecd_y0_isomorphic`

Distribution and root calculations:

- `ecd_solve`, `ecd_solve_sym`, `ecd_solve_trig`
- `ecd_solve_cusp_asym`, `ecd_y_slope`
- `ecd_pdf`, `ecd_cdf`, `ecd_ccdf`, `ecd_quantile`, `ecd_random`

Moments and transforms:

- `ecd_moment`, `ecd_statistics`, `ecd_asymptotic_statistics`
- `ecd_ellipticity`, `ecd_max_kurtosis`
- `ecd_cusp_std_moment`, `ecd_cusp_std_cf`, `ecd_cusp_std_mgf`
- `ecd_imgf`, `ecd_mu_d`, `ecd_ogf`, `ecd_rational`

## `ecld_models`

Type:

- `ecld_model`

Construction and distribution functions:

- `ecld_new`, `ecld_from_sd`, `ecld_solve`, `ecld_laplace_b`
- `ecld_const`, `ecld_pdf`, `ecld_cdf`, `ecld_ccdf`
- `ecld_moment`, `ecld_mean`, `ecld_variance`, `ecld_sd`
- `ecld_skewness`, `ecld_kurtosis`

Transforms and incomplete moments:

- `ecld_mgf`, `ecld_imgf`, `ecld_ogf`, `ecld_mu_d`
- `ecld_mgf_quartic`, `ecld_imgf_quartic`, `ecld_ogf_quartic`
- `ecld_mu_d_quartic`
- `ecld_incomplete_moment`, `ecld_imnt_sum`
- `ecld_ogf_star`, `ecld_ogf_star_hgeo`, `ecld_ogf_star_exp`
- `ecld_ogf_star_gamma_star`, `ecld_ogf_star_analytic`
- `ecld_gamma`, `ecld_gamma_hgeo`, `ecld_gamma_2f0`
- `ecld_y_slope`, `ecld_y_slope_trunc`

SGED wrappers:

- `ecld_sged_const`, `ecld_sged_cdf`, `ecld_sged_moment`
- `ecld_sged_mgf`, `ecld_sged_imgf`, `ecld_sged_ogf`

Option/lambda operators:

- `ecld_op_o`, `ecld_op_v`, `ecld_op_q`, `ecld_op_q_skew`
- `ecld_op_u_lag`, `ecld_op_vl_quartic`
- `ecld_quartic_q`, `ecld_quartic_qp`
- `ecld_quartic_qp_atm_ki`, `ecld_quartic_qp_rho`
- `ecld_quartic_qp_skew`, `ecld_quartic_qp_atm_skew`
- `ecld_quartic_sn0_atm_ki`, `ecld_quartic_sn0_rho_stdev`
- `ecld_quartic_sn0_skew`, `ecld_quartic_sn0_max_rnv`
- `ecld_quartic_model_sample`
- `ecld_fixed_point_atm_ki`, `ecld_fixed_point_shift`

## `ecd_processes`

Types:

- `sld_model`
- `sld_cumulant_result`

Laplace and standardized Lihn-Laplace:

- `laplace_pdf`, `laplace_random`
- `stdlap_pdf`, `stdlap_cdf`, `stdlap_quantile`, `stdlap_random`
- `stdlap_cf`, `stdlap_cumulants`, `stdlap_pdf_poly`

Stable count:

- `stable_pdf_positive`
- `stable_count_pdf`, `stable_count_cdf`, `stable_count_quantile`
- `stable_count_random`, `stable_count_cf`, `stable_count_cumulants`

SLD and QSLD:

- `sld_new`, `sld_pdf`, `sld_cdf`, `sld_quantile`, `sld_random`
- `sld_cf`, `sld_cumulants`
- `qsl_variance_analytic`, `qsl_skewness_analytic`
- `qsl_kurtosis_analytic`, `qsl_std_pdf0_analytic`
- `qsl_pdf_integrand_analytic`

Levy and moment utilities:

- `levy_dlambda`, `levy_dskewed`
- `k2moments`, `moments2k`

## `lamp_process`

Types:

- `lamp_model`
- `lamp_result`

Procedures:

- `lamp_new`, `lamp_sd_factor`, `stable_random`
- `lamp_generate_tau`, `lamp_stable_random_walk`
- `lamp_simulate_once`, `lamp_simulate`

## `ecd_options`

- `bs_option_price`, `bs_call_price`, `bs_put_price`
- `bs_implied_volatility`, `option_intrinsic_value`
- `polyfit_option`

## `ecd_fitting`

Types:

- `optimization_result`
- `ecd_fit_result`
- `ecld_fit_result`
- `sld_fit_result`

Procedures:

- `nelder_mead`
- `fit_ecd_mle`, compatibility name `ecd_standardfit`
- `fit_ecld_moments`, compatibility name `fit_ecld_mle`
- `fit_sld_mle`, compatibility name `qsld_fit`
- `ecd_estimate_const`

`fit_ecld_moments` retains its historical name, but it evaluates the full ECLD
likelihood after moment-based initialization.

## `ecd_timeseries`

Type:

- `sample_statistics`

Procedures:

- `difference_series`, `lag_series`
- `sample_stats`, `lag_stats`
- `empirical_quantile`, `quantilize`, `manage_hist_tails`

## `ecd_math`

- `normal_cdf`, `normal_quantile`
- `gamma_p`, `gamma_q`, `gamma_quantile`, `upper_incomplete_gamma`
- `integrate_adaptive`, `brent_root`
- `bessel_k`, `dawson_f`, `erfi_f`, `erfcx_f`
- `erfq`, `erfq_sum`
- `rational_approx`, `hypergeom_2f0`, `log_gamma_sign`

## `ecd_rng`

Type:

- `rng_state`

Procedures:

- `rng_seed`, `rng_uniform`, `rng_normal`
- `rng_exponential`, `rng_gamma`

## `ecd_compat`

Provides short names corresponding to important R exports, including `dec`,
`pec`, `qec`, `rec`, `dstdlap`, `pstdlap`, `qstdlap`, `rstdlap`, `dsl`, `psl`,
`qsl`, `rsl`, stable-count aliases, Black-Scholes aliases, and fitting aliases.
R names containing periods are represented with underscores or descriptive
Fortran names because periods are not valid in Fortran identifiers.
