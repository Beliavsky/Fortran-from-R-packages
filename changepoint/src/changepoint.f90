! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint
use r_kinds, only : dp
use changepoint_types, only : changepoint_result, amoc_result, binseg_result, segneigh_result, crops_solution
use changepoint_types, only : cp_ok, cp_invalid_argument, cp_invalid_data, cp_linalg_failure
use changepoint_costs, only : cp_cost_mean_normal, cp_cost_var_normal, cp_cost_meanvar_normal
use changepoint_costs, only : cp_cost_exponential, cp_cost_gamma, cp_cost_poisson
use changepoint_penalties, only : cp_penalty_value
use changepoint_decision, only : cp_decision
use changepoint_inference, only : cp_amoc_asymptotic_value
use changepoint_amoc, only : cp_amoc, cp_amoc_cusum, cp_amoc_css
use changepoint_multiple, only : cp_pelt, cp_binseg, cp_segneigh
use changepoint_nonparametric, only : cp_binseg_cusum, cp_binseg_css, cp_segneigh_cusum, cp_segneigh_css
use changepoint_crops, only : cp_crops
use changepoint_regression, only : cp_regression_amoc, cp_regression_pelt, cp_regression_segment_fit
use changepoint_fit, only : cp_segment_means, cp_segment_variances_mle, cp_segment_scales
use changepoint_fit, only : cp_segment_trend_fits
use changepoint_fit, only : cp_segment_regression_fits
implicit none
private

public :: dp
public :: changepoint_result, amoc_result, binseg_result, segneigh_result, crops_solution
public :: cp_ok, cp_invalid_argument, cp_invalid_data, cp_linalg_failure
public :: cp_cost_mean_normal, cp_cost_var_normal, cp_cost_meanvar_normal
public :: cp_cost_exponential, cp_cost_gamma, cp_cost_poisson
public :: cp_penalty_value
public :: cp_decision
public :: cp_amoc_asymptotic_value
public :: cp_amoc, cp_amoc_cusum, cp_amoc_css
public :: cp_pelt, cp_binseg, cp_segneigh
public :: cp_binseg_cusum, cp_binseg_css, cp_segneigh_cusum, cp_segneigh_css
public :: cp_crops
public :: cp_regression_amoc, cp_regression_pelt, cp_regression_segment_fit
public :: cp_segment_means, cp_segment_variances_mle, cp_segment_scales
public :: cp_segment_trend_fits
public :: cp_segment_regression_fits

end module changepoint
