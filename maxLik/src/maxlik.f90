! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik
  use maxlik_kinds, only: dp
  use maxlik_status
  use maxlik_types
  use maxlik_evaluation, only: initialize_problem, set_fixed, clear_fixed, set_bounds, &
    set_equality_constraints, set_inequality_constraints, clear_constraints, &
    numeric_gradient, numeric_hessian, constraint_violation
  use maxlik_solvers, only: solve_newton, solve_bfgs, solve_bfgsr, solve_bhhh, &
    solve_cg, solve_nelder_mead, solve_sann, solve_sga, solve_adam
  use maxlik_inference, only: covariance_matrix, robust_covariance_matrix, standard_errors, &
    normal_confidence_intervals, maxlik_aic, compare_derivatives, condition_number
  use maxlik_api, only: max_lik
  use maxlik_utilities
  implicit none
  public
end module maxlik
