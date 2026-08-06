! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif
  use rpeif_kinds, only : dp, pi
  use rpeif_types, only : nuisance_parameters, rpeif_options, influence_result, &
    rpeif_success, rpeif_invalid_argument, rpeif_numerical_failure, rpeif_unknown_estimator
  use rpeif_nuisance, only : nuisance_parameters_fn
  use rpeif_influence, only : influence_from_data, influence_from_nuisance, &
    influence_series, evaluate_shape, supported_estimator
  use rpeif_stats, only : lower_partial_moment, upper_partial_moment
  use rpeif_robust, only : robust_clean, robust_location_scale
  use rpeif_prewhiten, only : ar_prewhiten
  use rpeif_compat, only : if_mean, if_sd, if_semisd, if_var, if_es, if_sr, if_sortino, &
    if_downside_sharpe, if_es_ratio, if_var_ratio, if_rachev_ratio, if_robust_mean, &
    if_lpm, if_omega_ratio
  implicit none
  private
  public :: dp, pi
  public :: nuisance_parameters, rpeif_options, influence_result
  public :: rpeif_success, rpeif_invalid_argument, rpeif_numerical_failure, rpeif_unknown_estimator
  public :: nuisance_parameters_fn
  public :: influence_from_data, influence_from_nuisance, influence_series, evaluate_shape
  public :: supported_estimator
  public :: lower_partial_moment, upper_partial_moment
  public :: robust_clean, robust_location_scale, ar_prewhiten
  public :: if_mean, if_sd, if_semisd, if_var, if_es, if_sr, if_sortino
  public :: if_downside_sharpe, if_es_ratio, if_var_ratio, if_rachev_ratio
  public :: if_robust_mean, if_lpm, if_omega_ratio
end module rpeif
