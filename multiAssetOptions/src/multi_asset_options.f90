! Modern Fortran translation of the computational core of multiAssetOptions.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multi_asset_options
   use mao_kinds, only: dp
   use mao_status, only: status_type, mao_success, mao_invalid_argument, &
      mao_allocation_error, mao_solver_failure, mao_step_failure
   use mao_types, only: option_spec, fd_spec, time_spec, pricing_config, &
      asset_grid, grid_set, pricing_result, initialize_config, validate_config, &
      payoff_digital, payoff_best_of, payoff_worst_of, exercise_european, &
      exercise_american, option_call, option_put, timestep_constant, &
      timestep_adaptive
   use mao_sparse, only: csr_matrix, csr_matvec, csr_to_dense, csr_diagonal
   use mao_grid, only: node_spacer, build_grid, linear_index, decode_index, &
      interpolate_value
   use mao_payoff, only: payoff_values
   use mao_operator, only: build_fdm_operator
   use mao_pricing, only: price_multi_asset
   implicit none
   public

end module multi_asset_options
