! SPDX-License-Identifier: GPL-2.0-or-later
module highs_constants
   use highs_kinds, only : dp, highs_int
   implicit none
   private

   integer(highs_int), parameter, public :: highs_status_error = -1
   integer(highs_int), parameter, public :: highs_status_ok = 0
   integer(highs_int), parameter, public :: highs_status_warning = 1
   integer(highs_int), parameter, public :: highs_backend_unavailable = -100
   integer(highs_int), parameter, public :: highs_invalid_argument = -101

   integer(highs_int), parameter, public :: highs_minimize = 1
   integer(highs_int), parameter, public :: highs_maximize = -1

   integer(highs_int), parameter, public :: highs_matrix_colwise = 1
   integer(highs_int), parameter, public :: highs_matrix_rowwise = 2
   integer(highs_int), parameter, public :: highs_matrix_rowwise_partitioned = 3
   integer(highs_int), parameter, public :: highs_hessian_triangular = 1
   integer(highs_int), parameter, public :: highs_hessian_square = 2

   integer(highs_int), parameter, public :: highs_var_continuous = 0
   integer(highs_int), parameter, public :: highs_var_integer = 1
   integer(highs_int), parameter, public :: highs_var_semicontinuous = 2
   integer(highs_int), parameter, public :: highs_var_semiinteger = 3
   integer(highs_int), parameter, public :: highs_var_implicit_integer = 4

   integer(highs_int), parameter, public :: highs_basis_lower = 0
   integer(highs_int), parameter, public :: highs_basis_basic = 1
   integer(highs_int), parameter, public :: highs_basis_upper = 2
   integer(highs_int), parameter, public :: highs_basis_zero = 3
   integer(highs_int), parameter, public :: highs_basis_nonbasic = 4

   integer(highs_int), parameter, public :: highs_model_notset = 0
   integer(highs_int), parameter, public :: highs_model_load_error = 1
   integer(highs_int), parameter, public :: highs_model_error = 2
   integer(highs_int), parameter, public :: highs_model_presolve_error = 3
   integer(highs_int), parameter, public :: highs_model_solve_error = 4
   integer(highs_int), parameter, public :: highs_model_postsolve_error = 5
   integer(highs_int), parameter, public :: highs_model_empty = 6
   integer(highs_int), parameter, public :: highs_model_optimal = 7
   integer(highs_int), parameter, public :: highs_model_infeasible = 8
   integer(highs_int), parameter, public :: highs_model_unbounded_or_infeasible = 9
   integer(highs_int), parameter, public :: highs_model_unbounded = 10
   integer(highs_int), parameter, public :: highs_model_objective_bound = 11
   integer(highs_int), parameter, public :: highs_model_objective_target = 12
   integer(highs_int), parameter, public :: highs_model_time_limit = 13
   integer(highs_int), parameter, public :: highs_model_iteration_limit = 14
   integer(highs_int), parameter, public :: highs_model_unknown = 15
   integer(highs_int), parameter, public :: highs_model_solution_limit = 16
   integer(highs_int), parameter, public :: highs_model_interrupt = 17
   integer(highs_int), parameter, public :: highs_model_memory_limit = 18
   integer(highs_int), parameter, public :: highs_model_highs_interrupt = 19

   real(dp), parameter, public :: highs_default_infinity = 1.0e30_dp
end module highs_constants
