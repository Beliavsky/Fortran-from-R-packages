! SPDX-License-Identifier: GPL-2.0-or-later
module highs_solver_api
   use, intrinsic :: iso_c_binding, only : c_char, c_null_char, c_int, c_ptr, &
      c_null_ptr, c_associated
   use highs_kinds, only : dp, highs_int
   use highs_constants
   use highs_types
   use highs_c_bindings
   implicit none
   private

   public :: highs_backend_available, highs_load_backend, highs_unload_backend
   public :: highs_backend_error, highs_backend_version, highs_infinity
   public :: highs_new_solver, highs_destroy_solver, highs_pass_model
   public :: highs_set_control, highs_run, highs_get_solution
   public :: highs_set_bool_option, highs_set_int_option
   public :: highs_set_double_option, highs_set_string_option
   public :: highs_reset_options, highs_clear_model, highs_clear_solver
   public :: highs_presolve, highs_change_costs, highs_change_variable_bounds
   public :: highs_change_constraint_bounds, highs_change_coefficient
   public :: highs_change_integrality, highs_change_sense, highs_change_offset
   public :: highs_write_model, highs_read_model, highs_set_start
   public :: highs_get_basis, highs_set_basis, highs_clear_basis
   public :: highs_get_dual_ray, highs_get_primal_ray

contains

   logical function highs_backend_available() result(ok)
      ok = hf_api_backend_available() /= 0
   end function highs_backend_available

   subroutine highs_load_backend(path, status)
      character(len=*), intent(in), optional :: path
      integer(highs_int), intent(out), optional :: status
      character(c_char), allocatable :: cpath(:)
      integer(c_int) :: rc
      if (present(path)) then
         cpath = to_c_string(path)
      else
         allocate(cpath(1))
         cpath(1) = c_null_char
      end if
      rc = hf_api_load_backend(cpath)
      if (present(status)) status = merge(highs_status_ok, highs_backend_unavailable, rc /= 0)
   end subroutine highs_load_backend

   subroutine highs_unload_backend()
      call hf_api_unload_backend()
   end subroutine highs_unload_backend

   function highs_backend_error() result(message)
      character(len=:), allocatable :: message
      character(c_char) :: buffer(1024)
      integer(c_int) :: n
      n = hf_api_last_error(buffer, size(buffer, kind=c_int))
      message = from_c_buffer(buffer)
   end function highs_backend_error

   function highs_backend_version() result(version)
      character(len=:), allocatable :: version
      character(c_char) :: buffer(128)
      integer(c_int) :: rc
      rc = hf_api_version(buffer, size(buffer, kind=c_int))
      if (rc < 0) then
         version = "unavailable"
      else
         version = from_c_buffer(buffer)
      end if
   end function highs_backend_version

   real(dp) function highs_infinity(solver) result(value)
      type(highs_solver), intent(in), optional :: solver
      if (present(solver)) then
         value = hf_api_infinity(solver%handle)
      else
         value = hf_api_infinity(c_null_ptr)
      end if
   end function highs_infinity

   subroutine highs_new_solver(solver, status)
      type(highs_solver), intent(inout) :: solver
      integer(highs_int), intent(out), optional :: status
      integer(highs_int) :: stat
      if (c_associated(solver%handle)) call highs_destroy_solver(solver)
      if (.not. highs_backend_available()) then
         stat = highs_backend_unavailable
      else
         solver%handle = hf_api_create()
         if (c_associated(solver%handle)) then
            stat = highs_status_ok
         else
            stat = highs_status_error
         end if
      end if
      solver%num_col = 0
      solver%num_row = 0
      if (present(status)) status = stat
   end subroutine highs_new_solver

   subroutine highs_destroy_solver(solver)
      type(highs_solver), intent(inout) :: solver
      if (c_associated(solver%handle)) call hf_api_destroy(solver%handle)
      solver%handle = c_null_ptr
      solver%num_col = 0
      solver%num_row = 0
   end subroutine highs_destroy_solver

   subroutine highs_pass_model(solver, model, status)
      type(highs_solver), intent(inout) :: solver
      type(highs_model), intent(in) :: model
      integer(highs_int), intent(out) :: status
      integer(highs_int) :: rc
      if (.not. c_associated(solver%handle)) then
         call highs_new_solver(solver, rc)
         if (rc /= highs_status_ok) then
            status = rc
            return
         end if
      end if
      if (.not. model%valid()) then
         status = highs_invalid_argument
         return
      end if
      rc = hf_api_pass_model(solver%handle, model%num_col, model%num_row, &
         int(model%a%nnz(), highs_int), model%a%format, model%sense, model%offset, &
         model%col_cost, model%col_lower, model%col_upper, model%row_lower, &
         model%row_upper, model%a%start, model%a%index, model%a%value, model%integrality)
      if (rc == highs_status_error) then
         status = rc
         return
      end if
      if (model%has_hessian) then
         rc = hf_api_pass_hessian(solver%handle, model%num_col, &
            int(model%q%nnz(), highs_int), model%q%format, model%q%start, &
            model%q%index, model%q%value)
      end if
      solver%num_col = model%num_col
      solver%num_row = model%num_row
      status = rc
   end subroutine highs_pass_model

   subroutine highs_set_control(solver, control, status)
      type(highs_solver), intent(in) :: solver
      type(highs_control), intent(in) :: control
      integer(highs_int), intent(out) :: status
      integer(highs_int) :: rc
      status = highs_status_ok
      call highs_set_int_option(solver, "threads", control%threads, rc)
      call combine_status(status, rc)
      call highs_set_bool_option(solver, "log_to_console", control%log_to_console, rc)
      call combine_status(status, rc)
      call highs_set_bool_option(solver, "output_flag", control%log_to_console, rc)
      call combine_status(status, rc)
      call highs_set_string_option(solver, "solver", trim(control%solver), rc)
      call combine_status(status, rc)
      call highs_set_string_option(solver, "presolve", trim(control%presolve), rc)
      call combine_status(status, rc)
      call highs_set_string_option(solver, "parallel", &
         merge("on ", "off", control%parallel .or. control%threads > 1), rc)
      call combine_status(status, rc)
      if (control%time_limit >= 0.0_dp) then
         call highs_set_double_option(solver, "time_limit", control%time_limit, rc)
         call combine_status(status, rc)
      end if
      if (control%mip_rel_gap >= 0.0_dp) then
         call highs_set_double_option(solver, "mip_rel_gap", control%mip_rel_gap, rc)
         call combine_status(status, rc)
      end if
      if (control%mip_abs_gap >= 0.0_dp) then
         call highs_set_double_option(solver, "mip_abs_gap", control%mip_abs_gap, rc)
         call combine_status(status, rc)
      end if
      if (control%random_seed >= 0) then
         call highs_set_int_option(solver, "random_seed", control%random_seed, rc)
         call combine_status(status, rc)
      end if
   end subroutine highs_set_control

   subroutine highs_run(solver, status)
      type(highs_solver), intent(in) :: solver
      integer(highs_int), intent(out) :: status
      if (.not. c_associated(solver%handle)) then
         status = highs_invalid_argument
      else
         status = hf_api_run(solver%handle)
      end if
   end subroutine highs_run

   subroutine highs_get_solution(solver, solution, status)
      type(highs_solver), intent(in) :: solver
      type(highs_solution), intent(out) :: solution
      integer(highs_int), intent(out), optional :: status
      integer(c_int) :: rc, vv, dv
      type(hf_info_c) :: ci
      character(c_char) :: buffer(128)
      if (.not. c_associated(solver%handle)) then
         solution%call_status = highs_invalid_argument
         if (present(status)) status = solution%call_status
         return
      end if
      allocate(solution%col_value(solver%num_col), solution%col_dual(solver%num_col))
      allocate(solution%row_value(solver%num_row), solution%row_dual(solver%num_row))
      rc = hf_api_get_solution(solver%handle, solution%col_value, solution%col_dual, &
         solution%row_value, solution%row_dual, vv, dv)
      solution%call_status = rc
      solution%value_valid = vv /= 0
      solution%dual_valid = dv /= 0
      solution%model_status = hf_api_model_status(solver%handle)
      solution%objective_value = hf_api_objective_value(solver%handle)
      rc = hf_api_status_message(solver%handle, buffer, size(buffer, kind=c_int))
      solution%status_message = from_c_buffer(buffer)
      rc = hf_api_get_info(solver%handle, ci)
      call copy_info(ci, solution%info)
      if (present(status)) status = solution%call_status
   end subroutine highs_get_solution

   subroutine highs_set_bool_option(solver, name, value, status)
      type(highs_solver), intent(in) :: solver
      character(len=*), intent(in) :: name
      logical, intent(in) :: value
      integer(highs_int), intent(out) :: status
      character(c_char), allocatable :: cname(:)
      cname = to_c_string(name)
      status = hf_api_set_bool_option(solver%handle, cname, merge(1_c_int, 0_c_int, value))
   end subroutine highs_set_bool_option

   subroutine highs_set_int_option(solver, name, value, status)
      type(highs_solver), intent(in) :: solver
      character(len=*), intent(in) :: name
      integer(highs_int), intent(in) :: value
      integer(highs_int), intent(out) :: status
      character(c_char), allocatable :: cname(:)
      cname = to_c_string(name)
      status = hf_api_set_int_option(solver%handle, cname, value)
   end subroutine highs_set_int_option

   subroutine highs_set_double_option(solver, name, value, status)
      type(highs_solver), intent(in) :: solver
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      integer(highs_int), intent(out) :: status
      character(c_char), allocatable :: cname(:)
      cname = to_c_string(name)
      status = hf_api_set_double_option(solver%handle, cname, value)
   end subroutine highs_set_double_option

   subroutine highs_set_string_option(solver, name, value, status)
      type(highs_solver), intent(in) :: solver
      character(len=*), intent(in) :: name, value
      integer(highs_int), intent(out) :: status
      character(c_char), allocatable :: cname(:), cvalue(:)
      cname = to_c_string(name)
      cvalue = to_c_string(value)
      status = hf_api_set_string_option(solver%handle, cname, cvalue)
   end subroutine highs_set_string_option

   subroutine highs_reset_options(solver, status)
      type(highs_solver), intent(in) :: solver
      integer(highs_int), intent(out) :: status
      status = hf_api_reset_options(solver%handle)
   end subroutine highs_reset_options

   subroutine highs_clear_model(solver, status)
      type(highs_solver), intent(inout) :: solver
      integer(highs_int), intent(out) :: status
      status = hf_api_clear_model(solver%handle)
      if (status /= highs_status_error) then
         solver%num_col = 0
         solver%num_row = 0
      end if
   end subroutine highs_clear_model

   subroutine highs_clear_solver(solver, status)
      type(highs_solver), intent(in) :: solver
      integer(highs_int), intent(out) :: status
      status = hf_api_clear_solver(solver%handle)
   end subroutine highs_clear_solver

   subroutine highs_presolve(solver, status)
      type(highs_solver), intent(in) :: solver
      integer(highs_int), intent(out) :: status
      status = hf_api_presolve(solver%handle)
   end subroutine highs_presolve

   subroutine highs_change_costs(solver, index, value, status)
      type(highs_solver), intent(in) :: solver
      integer, intent(in) :: index(:)
      real(dp), intent(in) :: value(:)
      integer(highs_int), intent(out) :: status
      integer(highs_int), allocatable :: zero_index(:)
      if (size(index) /= size(value) .or. any(index < 1) .or. any(index > solver%num_col)) then
         status = highs_invalid_argument
         return
      end if
      zero_index = int(index - 1, highs_int)
      status = hf_api_change_cost(solver%handle, int(size(index), highs_int), zero_index, value)
   end subroutine highs_change_costs

   subroutine highs_change_variable_bounds(solver, index, lower, upper, status)
      type(highs_solver), intent(in) :: solver
      integer, intent(in) :: index(:)
      real(dp), intent(in) :: lower(:), upper(:)
      integer(highs_int), intent(out) :: status
      integer(highs_int), allocatable :: zero_index(:)
      if (size(index) /= size(lower) .or. size(index) /= size(upper) .or. &
          any(index < 1) .or. any(index > solver%num_col)) then
         status = highs_invalid_argument
         return
      end if
      zero_index = int(index - 1, highs_int)
      status = hf_api_change_col_bounds(solver%handle, int(size(index), highs_int), &
         zero_index, lower, upper)
   end subroutine highs_change_variable_bounds

   subroutine highs_change_constraint_bounds(solver, index, lower, upper, status)
      type(highs_solver), intent(in) :: solver
      integer, intent(in) :: index(:)
      real(dp), intent(in) :: lower(:), upper(:)
      integer(highs_int), intent(out) :: status
      integer(highs_int), allocatable :: zero_index(:)
      if (size(index) /= size(lower) .or. size(index) /= size(upper) .or. &
          any(index < 1) .or. any(index > solver%num_row)) then
         status = highs_invalid_argument
         return
      end if
      zero_index = int(index - 1, highs_int)
      status = hf_api_change_row_bounds(solver%handle, int(size(index), highs_int), &
         zero_index, lower, upper)
   end subroutine highs_change_constraint_bounds

   subroutine highs_change_coefficient(solver, row, col, value, status)
      type(highs_solver), intent(in) :: solver
      integer, intent(in) :: row, col
      real(dp), intent(in) :: value
      integer(highs_int), intent(out) :: status
      if (row < 1 .or. row > solver%num_row .or. col < 1 .or. col > solver%num_col) then
         status = highs_invalid_argument
      else
         status = hf_api_change_coeff(solver%handle, int(row - 1, highs_int), &
            int(col - 1, highs_int), value)
      end if
   end subroutine highs_change_coefficient

   subroutine highs_change_integrality(solver, index, vartype, status)
      type(highs_solver), intent(in) :: solver
      integer, intent(in) :: index(:)
      integer(highs_int), intent(in) :: vartype(:)
      integer(highs_int), intent(out) :: status
      integer(highs_int), allocatable :: zero_index(:)
      if (size(index) /= size(vartype) .or. any(index < 1) .or. &
          any(index > solver%num_col) .or. any(vartype < 0) .or. any(vartype > 4)) then
         status = highs_invalid_argument
         return
      end if
      zero_index = int(index - 1, highs_int)
      status = hf_api_change_integrality(solver%handle, int(size(index), highs_int), &
         zero_index, vartype)
   end subroutine highs_change_integrality

   subroutine highs_change_sense(solver, sense, status)
      type(highs_solver), intent(in) :: solver
      integer(highs_int), intent(in) :: sense
      integer(highs_int), intent(out) :: status
      if (sense /= highs_minimize .and. sense /= highs_maximize) then
         status = highs_invalid_argument
      else
         status = hf_api_change_sense(solver%handle, sense)
      end if
   end subroutine highs_change_sense

   subroutine highs_change_offset(solver, offset, status)
      type(highs_solver), intent(in) :: solver
      real(dp), intent(in) :: offset
      integer(highs_int), intent(out) :: status
      status = hf_api_change_offset(solver%handle, offset)
   end subroutine highs_change_offset

   subroutine highs_write_model(solver, filename, status)
      type(highs_solver), intent(in) :: solver
      character(len=*), intent(in) :: filename
      integer(highs_int), intent(out) :: status
      character(c_char), allocatable :: cfile(:)
      cfile = to_c_string(filename)
      status = hf_api_write_model(solver%handle, cfile)
   end subroutine highs_write_model

   subroutine highs_read_model(solver, filename, status)
      type(highs_solver), intent(inout) :: solver
      character(len=*), intent(in) :: filename
      integer(highs_int), intent(out) :: status
      character(c_char), allocatable :: cfile(:)
      if (.not. c_associated(solver%handle)) call highs_new_solver(solver, status)
      if (.not. c_associated(solver%handle)) return
      cfile = to_c_string(filename)
      status = hf_api_read_model(solver%handle, cfile)
      if (status /= highs_status_error) then
         solver%num_col = hf_api_num_col(solver%handle)
         solver%num_row = hf_api_num_row(solver%handle)
      end if
   end subroutine highs_read_model

   subroutine highs_set_start(solver, col_value, status, col_dual, row_value, row_dual)
      type(highs_solver), intent(in) :: solver
      real(dp), intent(in) :: col_value(:)
      integer(highs_int), intent(out) :: status
      real(dp), intent(in), optional :: col_dual(:), row_value(:), row_dual(:)
      real(dp), allocatable :: cd(:), rv(:), rd(:)
      integer(c_int) :: dual_valid
      if (size(col_value) /= solver%num_col) then
         status = highs_invalid_argument
         return
      end if
      allocate(cd(solver%num_col), source=0.0_dp)
      allocate(rv(solver%num_row), source=0.0_dp)
      allocate(rd(solver%num_row), source=0.0_dp)
      dual_valid = 0
      if (present(col_dual)) then
         if (size(col_dual) /= solver%num_col) then
            status = highs_invalid_argument
            return
         end if
         cd = col_dual
         dual_valid = 1
      end if
      if (present(row_value)) then
         if (size(row_value) /= solver%num_row) then
            status = highs_invalid_argument
            return
         end if
         rv = row_value
      end if
      if (present(row_dual)) then
         if (size(row_dual) /= solver%num_row) then
            status = highs_invalid_argument
            return
         end if
         rd = row_dual
         dual_valid = 1
      end if
      status = hf_api_set_solution(solver%handle, col_value, cd, rv, rd, 1_c_int, dual_valid)
   end subroutine highs_set_start

   subroutine highs_get_basis(solver, basis, status)
      type(highs_solver), intent(in) :: solver
      type(highs_basis), intent(out) :: basis
      integer(highs_int), intent(out) :: status
      integer(c_int) :: valid
      allocate(basis%col_status(solver%num_col), basis%row_status(solver%num_row))
      status = hf_api_get_basis(solver%handle, basis%col_status, basis%row_status, valid)
      basis%valid = valid /= 0
   end subroutine highs_get_basis

   subroutine highs_set_basis(solver, basis, status)
      type(highs_solver), intent(in) :: solver
      type(highs_basis), intent(in) :: basis
      integer(highs_int), intent(out) :: status
      if (.not. allocated(basis%col_status) .or. .not. allocated(basis%row_status)) then
         status = highs_invalid_argument
      else if (size(basis%col_status) /= solver%num_col .or. &
               size(basis%row_status) /= solver%num_row) then
         status = highs_invalid_argument
      else
         status = hf_api_set_basis(solver%handle, basis%col_status, basis%row_status)
      end if
   end subroutine highs_set_basis

   subroutine highs_clear_basis(solver, status)
      type(highs_solver), intent(in) :: solver
      integer(highs_int), intent(out) :: status
      status = hf_api_clear_basis(solver%handle)
   end subroutine highs_clear_basis

   subroutine highs_get_dual_ray(solver, ray, has_ray, status)
      type(highs_solver), intent(in) :: solver
      real(dp), allocatable, intent(out) :: ray(:)
      logical, intent(out) :: has_ray
      integer(highs_int), intent(out) :: status
      integer(c_int) :: h
      allocate(ray(solver%num_row), source=0.0_dp)
      status = hf_api_get_dual_ray(solver%handle, ray, h)
      has_ray = h /= 0
   end subroutine highs_get_dual_ray

   subroutine highs_get_primal_ray(solver, ray, has_ray, status)
      type(highs_solver), intent(in) :: solver
      real(dp), allocatable, intent(out) :: ray(:)
      logical, intent(out) :: has_ray
      integer(highs_int), intent(out) :: status
      integer(c_int) :: h
      allocate(ray(solver%num_col), source=0.0_dp)
      status = hf_api_get_primal_ray(solver%handle, ray, h)
      has_ray = h /= 0
   end subroutine highs_get_primal_ray

   function to_c_string(text) result(cstr)
      character(len=*), intent(in) :: text
      character(c_char), allocatable :: cstr(:)
      integer :: i, n
      n = len_trim(text)
      allocate(cstr(n + 1))
      do i = 1, n
         cstr(i) = text(i:i)
      end do
      cstr(n + 1) = c_null_char
   end function to_c_string

   function from_c_buffer(buffer) result(text)
      character(c_char), intent(in) :: buffer(:)
      character(len=:), allocatable :: text
      integer :: i, n
      n = 0
      do i = 1, size(buffer)
         if (buffer(i) == c_null_char) exit
         n = n + 1
      end do
      allocate(character(len=n) :: text)
      do i = 1, n
         text(i:i) = buffer(i)
      end do
   end function from_c_buffer

   subroutine copy_info(source, target)
      type(hf_info_c), intent(in) :: source
      type(highs_info), intent(out) :: target
      target%valid = source%valid /= 0
      target%mip_node_count = source%mip_node_count
      target%simplex_iteration_count = source%simplex_iteration_count
      target%ipm_iteration_count = source%ipm_iteration_count
      target%crossover_iteration_count = source%crossover_iteration_count
      target%qp_iteration_count = source%qp_iteration_count
      target%primal_solution_status = source%primal_solution_status
      target%dual_solution_status = source%dual_solution_status
      target%basis_validity = source%basis_validity
      target%objective_function_value = source%objective_function_value
      target%mip_dual_bound = source%mip_dual_bound
      target%mip_gap = source%mip_gap
      target%max_integrality_violation = source%max_integrality_violation
      target%num_primal_infeasibilities = source%num_primal_infeasibilities
      target%max_primal_infeasibility = source%max_primal_infeasibility
      target%sum_primal_infeasibilities = source%sum_primal_infeasibilities
      target%num_dual_infeasibilities = source%num_dual_infeasibilities
      target%max_dual_infeasibility = source%max_dual_infeasibility
      target%sum_dual_infeasibilities = source%sum_dual_infeasibilities
      target%run_time = source%run_time
   end subroutine copy_info

   subroutine combine_status(accumulator, value)
      integer(highs_int), intent(inout) :: accumulator
      integer(highs_int), intent(in) :: value
      if (value == highs_status_error .or. value == highs_backend_unavailable) then
         accumulator = value
      else if (value == highs_status_warning .and. accumulator == highs_status_ok) then
         accumulator = value
      end if
   end subroutine combine_status

end module highs_solver_api
