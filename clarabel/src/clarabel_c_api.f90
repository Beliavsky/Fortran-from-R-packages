! SPDX-License-Identifier: Apache-2.0
module clarabel_c_api
   use, intrinsic :: iso_c_binding
   implicit none
   private

   type, bind(c), public :: clarabel_csc_c
      integer(c_size_t) :: nrows
      integer(c_size_t) :: ncols
      integer(c_size_t) :: nnz
      type(c_ptr) :: colptr
      type(c_ptr) :: rowind
      type(c_ptr) :: values
   end type clarabel_csc_c

   type, bind(c), public :: clarabel_cone_c
      integer(c_int8_t) :: tag
      integer(c_size_t) :: dim
      real(c_double) :: parameter
      type(c_ptr) :: alpha
      integer(c_size_t) :: alpha_len
   end type clarabel_cone_c

   type, bind(c), public :: clarabel_settings_c
      integer(c_int32_t) :: max_iter
      real(c_double) :: time_limit
      integer(c_int32_t) :: verbose
      real(c_double) :: max_step_fraction
      real(c_double) :: tol_gap_abs, tol_gap_rel, tol_feas
      real(c_double) :: tol_infeas_abs, tol_infeas_rel, tol_ktratio
      real(c_double) :: reduced_tol_gap_abs, reduced_tol_gap_rel, reduced_tol_feas
      real(c_double) :: reduced_tol_infeas_abs, reduced_tol_infeas_rel, reduced_tol_ktratio
      integer(c_int32_t) :: equilibrate_enable, equilibrate_max_iter
      real(c_double) :: equilibrate_min_scaling, equilibrate_max_scaling
      real(c_double) :: linesearch_backtrack_step, min_switch_step_length, min_terminate_step_length
      integer(c_int32_t) :: max_threads, direct_kkt_solver, direct_solve_method
      integer(c_int32_t) :: static_regularization_enable
      real(c_double) :: static_regularization_constant, static_regularization_proportional
      integer(c_int32_t) :: dynamic_regularization_enable
      real(c_double) :: dynamic_regularization_eps, dynamic_regularization_delta
      integer(c_int32_t) :: iterative_refinement_enable
      real(c_double) :: iterative_refinement_reltol, iterative_refinement_abstol
      integer(c_int32_t) :: iterative_refinement_max_iter
      real(c_double) :: iterative_refinement_stop_ratio
      integer(c_int32_t) :: presolve_enable, input_sparse_dropzeros
      integer(c_int32_t) :: chordal_decomposition_enable, chordal_decomposition_merge_method
      integer(c_int32_t) :: chordal_decomposition_compact, chordal_decomposition_complete_dual
   end type clarabel_settings_c

   type, bind(c), public :: clarabel_result_c
      integer(c_int32_t) :: status, iterations
      real(c_double) :: obj_val, obj_val_dual, solve_time, r_prim, r_dual
      real(c_double) :: mu, sigma, step_length
      real(c_double) :: cost_primal, cost_dual, res_primal, res_dual
      real(c_double) :: res_primal_inf, res_dual_inf, gap_abs, gap_rel, ktratio
      integer(c_size_t) :: linear_solver_threads, linear_solver_nnz_a, linear_solver_nnz_l
   end type clarabel_result_c

   public :: c_clarabel_settings_default, c_clarabel_solver_create
   public :: c_clarabel_solver_solve, c_clarabel_solver_update
   public :: c_clarabel_solver_is_update_allowed, c_clarabel_solver_free

   interface
      subroutine c_clarabel_settings_default(settings) bind(c, name="clarabel_settings_default")
         import :: clarabel_settings_c
         type(clarabel_settings_c), intent(out) :: settings
      end subroutine c_clarabel_settings_default

      function c_clarabel_solver_create(p, q, q_len, a, b, b_len, cones, ncones, settings, &
                                        out, error, error_capacity) result(code) &
                                        bind(c, name="clarabel_solver_create")
         import :: c_int32_t, c_size_t, c_ptr, c_char, clarabel_csc_c, clarabel_cone_c, clarabel_settings_c
         type(clarabel_csc_c), intent(in) :: p, a
         type(c_ptr), value :: q, b, cones
         integer(c_size_t), value :: q_len, b_len, ncones
         type(clarabel_settings_c), intent(in) :: settings
         type(c_ptr), intent(out) :: out
         character(c_char), intent(out) :: error(*)
         integer(c_size_t), value :: error_capacity
         integer(c_int32_t) :: code
      end function c_clarabel_solver_create

      function c_clarabel_solver_solve(solver, x, x_len, z, z_len, s, s_len, result, error, &
                                       error_capacity) result(code) bind(c, name="clarabel_solver_solve")
         import :: c_int32_t, c_size_t, c_ptr, c_char, clarabel_result_c
         type(c_ptr), value :: solver, x, z, s
         integer(c_size_t), value :: x_len, z_len, s_len
         type(clarabel_result_c), intent(out) :: result
         character(c_char), intent(out) :: error(*)
         integer(c_size_t), value :: error_capacity
         integer(c_int32_t) :: code
      end function c_clarabel_solver_solve

      function c_clarabel_solver_update(solver, p_values, p_len, a_values, a_len, q, q_len, b, b_len, &
                                        error, error_capacity) result(code) bind(c, name="clarabel_solver_update")
         import :: c_int32_t, c_size_t, c_ptr, c_char
         type(c_ptr), value :: solver, p_values, a_values, q, b
         integer(c_size_t), value :: p_len, a_len, q_len, b_len
         character(c_char), intent(out) :: error(*)
         integer(c_size_t), value :: error_capacity
         integer(c_int32_t) :: code
      end function c_clarabel_solver_update

      function c_clarabel_solver_is_update_allowed(solver) result(flag) &
            bind(c, name="clarabel_solver_is_update_allowed")
         import :: c_int32_t, c_ptr
         type(c_ptr), value :: solver
         integer(c_int32_t) :: flag
      end function c_clarabel_solver_is_update_allowed

      subroutine c_clarabel_solver_free(solver) bind(c, name="clarabel_solver_free")
         import :: c_ptr
         type(c_ptr), value :: solver
      end subroutine c_clarabel_solver_free
   end interface

end module clarabel_c_api
