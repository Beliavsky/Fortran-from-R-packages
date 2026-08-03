! SPDX-License-Identifier: Apache-2.0
module osqp_c_bindings
   use, intrinsic :: iso_c_binding
   implicit none
   private

   type, bind(c), public :: of_settings_c
      integer(c_int) :: device
      integer(c_int) :: linsys_solver
      integer(c_int) :: allocate_solution
      integer(c_int) :: verbose
      integer(c_int) :: profiler_level
      integer(c_int) :: warm_starting
      integer(c_int) :: scaling
      integer(c_int) :: polishing
      real(c_double) :: rho
      integer(c_int) :: rho_is_vec
      real(c_double) :: sigma
      real(c_double) :: alpha
      integer(c_int) :: cg_max_iter
      integer(c_int) :: cg_tol_reduction
      real(c_double) :: cg_tol_fraction
      integer(c_int) :: cg_precond
      integer(c_int) :: adaptive_rho
      integer(c_int) :: adaptive_rho_interval
      real(c_double) :: adaptive_rho_fraction
      real(c_double) :: adaptive_rho_tolerance
      integer(c_int) :: max_iter
      real(c_double) :: eps_abs
      real(c_double) :: eps_rel
      real(c_double) :: eps_prim_inf
      real(c_double) :: eps_dual_inf
      integer(c_int) :: scaled_termination
      integer(c_int) :: check_termination
      integer(c_int) :: check_dualgap
      real(c_double) :: time_limit
      real(c_double) :: delta
      integer(c_int) :: polish_refine_iter
   end type of_settings_c

   type, bind(c), public :: of_info_c
      integer(c_int) :: status_val
      integer(c_int) :: status_polish
      real(c_double) :: obj_val
      real(c_double) :: dual_obj_val
      real(c_double) :: prim_res
      real(c_double) :: dual_res
      real(c_double) :: duality_gap
      integer(c_int) :: iter
      integer(c_int) :: rho_updates
      real(c_double) :: rho_estimate
      real(c_double) :: setup_time
      real(c_double) :: solve_time
      real(c_double) :: update_time
      real(c_double) :: polish_time
      real(c_double) :: run_time
      real(c_double) :: primdual_int
      real(c_double) :: rel_kkt_error
   end type of_info_c

   public :: c_of_api_load_backend, c_of_api_backend_available, c_of_api_unload_backend
   public :: c_of_api_last_error, c_of_api_create, c_of_api_destroy, c_of_api_solve
   public :: c_of_api_get_dimensions, c_of_api_get_solution, c_of_api_status_message
   public :: c_of_api_version, c_of_api_capabilities, c_of_api_infinity
   public :: c_of_api_error_message, c_of_api_get_settings
   public :: c_of_api_update_data_vec, c_of_api_update_data_mat
   public :: c_of_api_warm_start, c_of_api_cold_start, c_of_api_update_settings

   interface
      function c_of_api_load_backend(path) bind(c, name="of_api_load_backend") result(ok)
         import :: c_ptr, c_int
         type(c_ptr), value :: path
         integer(c_int) :: ok
      end function c_of_api_load_backend

      function c_of_api_backend_available() bind(c, name="of_api_backend_available") result(ok)
         import :: c_int
         integer(c_int) :: ok
      end function c_of_api_backend_available

      subroutine c_of_api_unload_backend() bind(c, name="of_api_unload_backend")
      end subroutine c_of_api_unload_backend

      function c_of_api_last_error(buffer, n) bind(c, name="of_api_last_error") result(length)
         import :: c_ptr, c_int
         type(c_ptr), value :: buffer
         integer(c_int), value :: n
         integer(c_int) :: length
      end function c_of_api_last_error

      function c_of_api_create(n, m, p_nnz, p_colptr, p_rowind, p_value, &
            a_nnz, a_colptr, a_rowind, a_value, q, l, u, settings, error_code) &
            bind(c, name="of_api_create") result(handle)
         import :: c_int, c_ptr, of_settings_c
         integer(c_int), value :: n, m, p_nnz, a_nnz
         type(c_ptr), value :: p_colptr, p_rowind, p_value
         type(c_ptr), value :: a_colptr, a_rowind, a_value
         type(c_ptr), value :: q, l, u
         type(of_settings_c), intent(in) :: settings
         integer(c_int), intent(out) :: error_code
         type(c_ptr) :: handle
      end function c_of_api_create

      subroutine c_of_api_destroy(handle) bind(c, name="of_api_destroy")
         import :: c_ptr
         type(c_ptr), value :: handle
      end subroutine c_of_api_destroy

      function c_of_api_solve(handle) bind(c, name="of_api_solve") result(status)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int) :: status
      end function c_of_api_solve

      function c_of_api_get_dimensions(handle, n, m) bind(c, name="of_api_get_dimensions") result(status)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), intent(out) :: n, m
         integer(c_int) :: status
      end function c_of_api_get_dimensions

      function c_of_api_get_solution(handle, x, y, prim_cert, dual_cert, info) &
            bind(c, name="of_api_get_solution") result(status)
         import :: c_ptr, c_int, of_info_c
         type(c_ptr), value :: handle, x, y, prim_cert, dual_cert
         type(of_info_c), intent(out) :: info
         integer(c_int) :: status
      end function c_of_api_get_solution

      function c_of_api_status_message(handle, buffer, n) bind(c, name="of_api_status_message") result(length)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle, buffer
         integer(c_int), value :: n
         integer(c_int) :: length
      end function c_of_api_status_message

      function c_of_api_version(buffer, n) bind(c, name="of_api_version") result(length)
         import :: c_ptr, c_int
         type(c_ptr), value :: buffer
         integer(c_int), value :: n
         integer(c_int) :: length
      end function c_of_api_version

      function c_of_api_capabilities() bind(c, name="of_api_capabilities") result(value)
         import :: c_int
         integer(c_int) :: value
      end function c_of_api_capabilities

      function c_of_api_infinity() bind(c, name="of_api_infinity") result(value)
         import :: c_double
         real(c_double) :: value
      end function c_of_api_infinity

      function c_of_api_error_message(code, buffer, n) bind(c, name="of_api_error_message") result(length)
         import :: c_ptr, c_int
         integer(c_int), value :: code, n
         type(c_ptr), value :: buffer
         integer(c_int) :: length
      end function c_of_api_error_message

      function c_of_api_get_settings(handle, settings) bind(c, name="of_api_get_settings") result(status)
         import :: c_ptr, c_int, of_settings_c
         type(c_ptr), value :: handle
         type(of_settings_c), intent(out) :: settings
         integer(c_int) :: status
      end function c_of_api_get_settings

      function c_of_api_update_data_vec(handle, q, has_q, l, has_l, u, has_u) &
            bind(c, name="of_api_update_data_vec") result(status)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle, q, l, u
         integer(c_int), value :: has_q, has_l, has_u
         integer(c_int) :: status
      end function c_of_api_update_data_vec

      function c_of_api_update_data_mat(handle, px, pidx, pn, has_p, ax, aidx, an, has_a) &
            bind(c, name="of_api_update_data_mat") result(status)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle, px, pidx, ax, aidx
         integer(c_int), value :: pn, has_p, an, has_a
         integer(c_int) :: status
      end function c_of_api_update_data_mat

      function c_of_api_warm_start(handle, x, has_x, y, has_y) &
            bind(c, name="of_api_warm_start") result(status)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle, x, y
         integer(c_int), value :: has_x, has_y
         integer(c_int) :: status
      end function c_of_api_warm_start

      function c_of_api_cold_start(handle) bind(c, name="of_api_cold_start") result(status)
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int) :: status
      end function c_of_api_cold_start

      function c_of_api_update_settings(handle, settings) bind(c, name="of_api_update_settings") result(status)
         import :: c_ptr, c_int, of_settings_c
         type(c_ptr), value :: handle
         type(of_settings_c), intent(in) :: settings
         integer(c_int) :: status
      end function c_of_api_update_settings
   end interface

end module osqp_c_bindings
