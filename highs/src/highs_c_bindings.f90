! SPDX-License-Identifier: GPL-2.0-or-later
module highs_c_bindings
   use, intrinsic :: iso_c_binding
   implicit none
   private

   type, bind(c), public :: hf_info_c
      integer(c_int) :: valid
      integer(c_int64_t) :: mip_node_count
      integer(c_int) :: simplex_iteration_count
      integer(c_int) :: ipm_iteration_count
      integer(c_int) :: crossover_iteration_count
      integer(c_int) :: qp_iteration_count
      integer(c_int) :: primal_solution_status
      integer(c_int) :: dual_solution_status
      integer(c_int) :: basis_validity
      real(c_double) :: objective_function_value
      real(c_double) :: mip_dual_bound
      real(c_double) :: mip_gap
      real(c_double) :: max_integrality_violation
      integer(c_int) :: num_primal_infeasibilities
      real(c_double) :: max_primal_infeasibility
      real(c_double) :: sum_primal_infeasibilities
      integer(c_int) :: num_dual_infeasibilities
      real(c_double) :: max_dual_infeasibility
      real(c_double) :: sum_dual_infeasibilities
      real(c_double) :: run_time
   end type hf_info_c

   public :: hf_api_load_backend, hf_api_backend_available, hf_api_unload_backend
   public :: hf_api_last_error, hf_api_create, hf_api_destroy, hf_api_pass_model
   public :: hf_api_pass_hessian, hf_api_run, hf_api_model_status
   public :: hf_api_objective_value, hf_api_get_solution, hf_api_get_info
   public :: hf_api_status_message, hf_api_version, hf_api_infinity
   public :: hf_api_num_col, hf_api_num_row, hf_api_get_sense
   public :: hf_api_set_bool_option, hf_api_set_int_option
   public :: hf_api_set_double_option, hf_api_set_string_option
   public :: hf_api_reset_options, hf_api_clear_model, hf_api_clear_solver
   public :: hf_api_presolve, hf_api_change_cost, hf_api_change_col_bounds
   public :: hf_api_change_row_bounds, hf_api_change_coeff
   public :: hf_api_change_integrality, hf_api_change_sense, hf_api_change_offset
   public :: hf_api_write_model, hf_api_read_model, hf_api_set_solution
   public :: hf_api_get_basis, hf_api_set_basis, hf_api_clear_basis
   public :: hf_api_get_dual_ray, hf_api_get_primal_ray

   interface
      integer(c_int) function hf_api_load_backend(path) bind(c)
         import c_int, c_char
         character(c_char), intent(in) :: path(*)
      end function
      integer(c_int) function hf_api_backend_available() bind(c)
         import c_int
      end function
      subroutine hf_api_unload_backend() bind(c)
      end subroutine
      integer(c_int) function hf_api_last_error(buffer, n) bind(c)
         import c_int, c_char
         character(c_char), intent(out) :: buffer(*)
         integer(c_int), value :: n
      end function
      type(c_ptr) function hf_api_create() bind(c)
         import c_ptr
      end function
      subroutine hf_api_destroy(handle) bind(c)
         import c_ptr
         type(c_ptr), value :: handle
      end subroutine
      integer(c_int) function hf_api_pass_model(handle, nc, nr, nz, fmt, sense, offset, &
            cost, cl, cu, rl, ru, start, index, value, integrality) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: nc, nr, nz, fmt, sense
         real(c_double), value :: offset
         real(c_double), intent(in) :: cost(*), cl(*), cu(*), rl(*), ru(*), value(*)
         integer(c_int), intent(in) :: start(*), index(*), integrality(*)
      end function
      integer(c_int) function hf_api_pass_hessian(handle, n, nz, fmt, start, index, value) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n, nz, fmt
         integer(c_int), intent(in) :: start(*), index(*)
         real(c_double), intent(in) :: value(*)
      end function
      integer(c_int) function hf_api_run(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_model_status(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      real(c_double) function hf_api_objective_value(handle) bind(c)
         import c_ptr, c_double
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_get_solution(handle, cv, cd, rv, rd, vv, dv) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: cv(*), cd(*), rv(*), rd(*)
         integer(c_int), intent(out) :: vv, dv
      end function
      integer(c_int) function hf_api_get_info(handle, info) bind(c)
         import c_ptr, c_int, hf_info_c
         type(c_ptr), value :: handle
         type(hf_info_c), intent(out) :: info
      end function
      integer(c_int) function hf_api_status_message(handle, buffer, n) bind(c)
         import c_ptr, c_int, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(out) :: buffer(*)
         integer(c_int), value :: n
      end function
      integer(c_int) function hf_api_version(buffer, n) bind(c)
         import c_int, c_char
         character(c_char), intent(out) :: buffer(*)
         integer(c_int), value :: n
      end function
      real(c_double) function hf_api_infinity(handle) bind(c)
         import c_ptr, c_double
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_num_col(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_num_row(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_get_sense(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_set_bool_option(handle, name, value) bind(c)
         import c_ptr, c_int, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: name(*)
         integer(c_int), value :: value
      end function
      integer(c_int) function hf_api_set_int_option(handle, name, value) bind(c)
         import c_ptr, c_int, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: name(*)
         integer(c_int), value :: value
      end function
      integer(c_int) function hf_api_set_double_option(handle, name, value) bind(c)
         import c_ptr, c_int, c_double, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: name(*)
         real(c_double), value :: value
      end function
      integer(c_int) function hf_api_set_string_option(handle, name, value) bind(c)
         import c_ptr, c_int, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: name(*), value(*)
      end function
      integer(c_int) function hf_api_reset_options(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_clear_model(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_clear_solver(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_presolve(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_change_cost(handle, n, index, value) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         integer(c_int), intent(in) :: index(*)
         real(c_double), intent(in) :: value(*)
      end function
      integer(c_int) function hf_api_change_col_bounds(handle, n, index, lower, upper) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         integer(c_int), intent(in) :: index(*)
         real(c_double), intent(in) :: lower(*), upper(*)
      end function
      integer(c_int) function hf_api_change_row_bounds(handle, n, index, lower, upper) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         integer(c_int), intent(in) :: index(*)
         real(c_double), intent(in) :: lower(*), upper(*)
      end function
      integer(c_int) function hf_api_change_coeff(handle, row, col, value) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: row, col
         real(c_double), value :: value
      end function
      integer(c_int) function hf_api_change_integrality(handle, n, index, value) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         integer(c_int), intent(in) :: index(*), value(*)
      end function
      integer(c_int) function hf_api_change_sense(handle, sense) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), value :: sense
      end function
      integer(c_int) function hf_api_change_offset(handle, value) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         real(c_double), value :: value
      end function
      integer(c_int) function hf_api_write_model(handle, filename) bind(c)
         import c_ptr, c_int, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: filename(*)
      end function
      integer(c_int) function hf_api_read_model(handle, filename) bind(c)
         import c_ptr, c_int, c_char
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: filename(*)
      end function
      integer(c_int) function hf_api_set_solution(handle, cv, cd, rv, rd, vv, dv) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         real(c_double), intent(in) :: cv(*), cd(*), rv(*), rd(*)
         integer(c_int), value :: vv, dv
      end function
      integer(c_int) function hf_api_get_basis(handle, cs, rs, valid) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: cs(*), rs(*), valid
      end function
      integer(c_int) function hf_api_set_basis(handle, cs, rs) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), intent(in) :: cs(*), rs(*)
      end function
      integer(c_int) function hf_api_clear_basis(handle) bind(c)
         import c_ptr, c_int
         type(c_ptr), value :: handle
      end function
      integer(c_int) function hf_api_get_dual_ray(handle, value, has_ray) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: value(*)
         integer(c_int), intent(out) :: has_ray
      end function
      integer(c_int) function hf_api_get_primal_ray(handle, value, has_ray) bind(c)
         import c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         real(c_double), intent(out) :: value(*)
         integer(c_int), intent(out) :: has_ray
      end function
   end interface
end module highs_c_bindings
