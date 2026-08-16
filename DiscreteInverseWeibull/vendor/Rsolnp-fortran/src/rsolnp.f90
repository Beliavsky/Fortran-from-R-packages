! SPDX-License-Identifier: GPL-2.0-only
module rsolnp
   use rsolnp_kinds, only : dp
   use rsolnp_callbacks, only : objective_callback, gradient_callback, vector_callback, jacobian_callback
   use rsolnp_types, only : solnp_control, solnp_problem, solnp_result, multistart_result, &
      kkt_diagnostics, problem_table_entry, solnp_success, solnp_max_iterations, &
      solnp_invalid_problem, solnp_numerical_failure, solnp_infeasible
   use rsolnp_problem, only : prepare_problem, solnp_standardize_problem
   use rsolnp_solver, only : solnp, csolnp, kkt_diagnose
   use rsolnp_multistart, only : startpars, csolnp_ms, gosolnp
   use rsolnp_benchmarks, only : solnp_problem_suite, solnp_problems_table
   implicit none
   private

   public :: dp
   public :: objective_callback, gradient_callback, vector_callback, jacobian_callback
   public :: solnp_control, solnp_problem, solnp_result, multistart_result
   public :: kkt_diagnostics, problem_table_entry
   public :: solnp_success, solnp_max_iterations, solnp_invalid_problem
   public :: solnp_numerical_failure, solnp_infeasible
   public :: prepare_problem, solnp_standardize_problem
   public :: solnp, csolnp, kkt_diagnose
   public :: startpars, csolnp_ms, gosolnp
   public :: solnp_problem_suite, solnp_problems_table

end module rsolnp
