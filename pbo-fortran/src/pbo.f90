! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo
  use pbo_kinds, only : dp
  use pbo_types, only : pbo_result, dominance_result, selection_result
  use pbo_core, only : performance_function, compute_pbo
  use pbo_metrics, only : column_mean, column_sum, sharpe_ratio, sharpe_ratio_rf, &
    omega_ratio, omega_ratio_threshold
  use pbo_analysis, only : dominance_curve, selection_frequencies
  use pbo_combinations, only : binomial_coefficient, generate_combinations
  implicit none
  public
end module pbo
