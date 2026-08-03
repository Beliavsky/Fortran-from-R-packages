! SPDX-License-Identifier: Apache-2.0
module osqp_constants
   use osqp_kinds, only : osqp_int, dp
   implicit none
   private

   integer(osqp_int), parameter, public :: osqp_backend_unavailable = -100
   integer(osqp_int), parameter, public :: osqp_invalid_argument = -101

   integer(osqp_int), parameter, public :: osqp_no_error = 0
   integer(osqp_int), parameter, public :: osqp_data_validation_error = 1
   integer(osqp_int), parameter, public :: osqp_settings_validation_error = 2
   integer(osqp_int), parameter, public :: osqp_linsys_solver_init_error = 3
   integer(osqp_int), parameter, public :: osqp_nonconvex_error = 4
   integer(osqp_int), parameter, public :: osqp_mem_alloc_error = 5
   integer(osqp_int), parameter, public :: osqp_workspace_not_init_error = 6
   integer(osqp_int), parameter, public :: osqp_algebra_load_error = 7
   integer(osqp_int), parameter, public :: osqp_fopen_error = 8
   integer(osqp_int), parameter, public :: osqp_codegen_defines_error = 9
   integer(osqp_int), parameter, public :: osqp_data_not_initialized = 10
   integer(osqp_int), parameter, public :: osqp_func_not_implemented = 11

   integer(osqp_int), parameter, public :: osqp_solved = 1
   integer(osqp_int), parameter, public :: osqp_solved_inaccurate = 2
   integer(osqp_int), parameter, public :: osqp_primal_infeasible = 3
   integer(osqp_int), parameter, public :: osqp_primal_infeasible_inaccurate = 4
   integer(osqp_int), parameter, public :: osqp_dual_infeasible = 5
   integer(osqp_int), parameter, public :: osqp_dual_infeasible_inaccurate = 6
   integer(osqp_int), parameter, public :: osqp_max_iter_reached = 7
   integer(osqp_int), parameter, public :: osqp_time_limit_reached = 8
   integer(osqp_int), parameter, public :: osqp_non_cvx = 9
   integer(osqp_int), parameter, public :: osqp_sigint = 10
   integer(osqp_int), parameter, public :: osqp_unsolved = 11

   integer(osqp_int), parameter, public :: osqp_unknown_solver = 0
   integer(osqp_int), parameter, public :: osqp_direct_solver = 1
   integer(osqp_int), parameter, public :: osqp_indirect_solver = 2
   integer(osqp_int), parameter, public :: osqp_no_preconditioner = 0
   integer(osqp_int), parameter, public :: osqp_diagonal_preconditioner = 1

   integer(osqp_int), parameter, public :: osqp_adaptive_rho_disabled = 0
   integer(osqp_int), parameter, public :: osqp_adaptive_rho_iterations = 1
   integer(osqp_int), parameter, public :: osqp_adaptive_rho_time = 2
   integer(osqp_int), parameter, public :: osqp_adaptive_rho_kkt_error = 3

   integer(osqp_int), parameter, public :: osqp_capability_direct_solver = int(z'01', osqp_int)
   integer(osqp_int), parameter, public :: osqp_capability_indirect_solver = int(z'02', osqp_int)
   integer(osqp_int), parameter, public :: osqp_capability_codegen = int(z'04', osqp_int)
   integer(osqp_int), parameter, public :: osqp_capability_update_matrices = int(z'08', osqp_int)
   integer(osqp_int), parameter, public :: osqp_capability_derivatives = int(z'10', osqp_int)

   real(dp), parameter, public :: osqp_infinity_default = 1.0e30_dp
end module osqp_constants
