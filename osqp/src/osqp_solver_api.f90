! SPDX-License-Identifier: Apache-2.0
module osqp_solver_api
   use, intrinsic :: iso_c_binding
   use osqp_kinds, only : dp, osqp_int
   use osqp_constants
   use osqp_c_bindings
   use osqp_types
   implicit none
   private

   public :: osqp_load_backend, osqp_backend_available, osqp_unload_backend
   public :: osqp_backend_error, osqp_backend_version, osqp_capabilities
   public :: osqp_infinity, osqp_error_string
   public :: osqp_setup, osqp_cleanup, osqp_solve_solver
   public :: osqp_update, osqp_warm_start, osqp_cold_start
   public :: osqp_update_settings, osqp_get_settings, osqp_get_dimensions

contains

   subroutine osqp_load_backend(path, status)
      character(len=*), intent(in), optional :: path
      integer(osqp_int), intent(out), optional :: status
      character(kind=c_char), allocatable, target :: cpath(:)
      integer(c_int) :: rc
      integer :: i, n
      if (present(path)) then
         n = len_trim(path)
         allocate(cpath(n+1))
         do i = 1, n
            cpath(i) = char(iachar(path(i:i)), kind=c_char)
         end do
         cpath(n+1) = c_null_char
         rc = c_of_api_load_backend(c_loc(cpath(1)))
      else
         rc = c_of_api_load_backend(c_null_ptr)
      end if
      if (present(status)) status = merge(osqp_no_error, osqp_backend_unavailable, rc /= 0)
   end subroutine osqp_load_backend

   logical function osqp_backend_available() result(ok)
      ok = c_of_api_backend_available() /= 0
   end function osqp_backend_available

   subroutine osqp_unload_backend()
      call c_of_api_unload_backend()
   end subroutine osqp_unload_backend

   function osqp_backend_error() result(message)
      character(len=:), allocatable :: message
      character(kind=c_char), target :: buffer(1024)
      integer(c_int) :: n
      n = c_of_api_last_error(c_loc(buffer(1)), int(size(buffer), c_int))
      message = c_buffer_to_string(buffer, n)
   end function osqp_backend_error

   function osqp_backend_version() result(version)
      character(len=:), allocatable :: version
      character(kind=c_char), target :: buffer(128)
      integer(c_int) :: n
      if (.not. osqp_backend_available()) then
         version = "unavailable"
         return
      end if
      n = c_of_api_version(c_loc(buffer(1)), int(size(buffer), c_int))
      version = c_buffer_to_string(buffer, n)
   end function osqp_backend_version

   integer(osqp_int) function osqp_capabilities() result(value)
      if (osqp_backend_available()) then
         value = c_of_api_capabilities()
      else
         value = 0
      end if
   end function osqp_capabilities

   real(dp) function osqp_infinity() result(value)
      if (osqp_backend_available()) then
         value = c_of_api_infinity()
      else
         value = osqp_infinity_default
      end if
   end function osqp_infinity

   function osqp_error_string(code) result(message)
      integer(osqp_int), intent(in) :: code
      character(len=:), allocatable :: message
      character(kind=c_char), target :: buffer(256)
      integer(c_int) :: n
      if (code == osqp_backend_unavailable) then
         message = "OSQP backend unavailable: " // osqp_backend_error()
      else if (code == osqp_invalid_argument) then
         message = "invalid Fortran API argument"
      else if (.not. osqp_backend_available()) then
         message = "OSQP backend unavailable"
      else
         n = c_of_api_error_message(code, c_loc(buffer(1)), int(size(buffer), c_int))
         message = c_buffer_to_string(buffer, n)
      end if
   end function osqp_error_string

   subroutine osqp_setup(solver, model, status, settings)
      type(osqp_solver), intent(inout) :: solver
      type(osqp_model), intent(in) :: model
      integer(osqp_int), intent(out) :: status
      type(osqp_settings), intent(in), optional :: settings
      type(osqp_settings) :: cfg
      type(of_settings_c) :: cset
      integer(c_int), allocatable, target :: pp(:), pi(:), ap(:), ai(:)
      real(c_double), allocatable, target :: px(:), ax(:), q(:), l(:), u(:)
      type(c_ptr) :: ppi, pii, ppx, api, aii, apx, qp, lp, up
      integer(c_int) :: rc

      status = osqp_invalid_argument
      if (.not. model%valid()) return
      cfg = osqp_settings()
      if (present(settings)) cfg = settings
      if (.not. cfg%valid()) return
      if (.not. osqp_backend_available()) then
         status = osqp_backend_unavailable
         return
      end if
      call osqp_cleanup(solver)
      call model_to_c_arrays(model, pp, pi, px, ap, ai, ax, q, l, u)
      ppi = c_loc(pp(1))
      pii = pointer_or_null_int(pi)
      ppx = pointer_or_null_real(px)
      api = c_loc(ap(1))
      aii = pointer_or_null_int(ai)
      apx = pointer_or_null_real(ax)
      qp = c_loc(q(1))
      lp = pointer_or_null_real(l)
      up = pointer_or_null_real(u)
      cset = settings_to_c(cfg)
      rc = osqp_backend_unavailable
      solver%handle = c_of_api_create(model%n, model%m, size(px), ppi, pii, ppx, &
         size(ax), api, aii, apx, qp, lp, up, cset, rc)
      if (.not. c_associated(solver%handle)) then
         status = rc
         solver%initialized = .false.
         return
      end if
      solver%model = model
      solver%settings = cfg
      solver%initialized = .true.
      status = osqp_no_error
   end subroutine osqp_setup

   subroutine osqp_cleanup(solver)
      type(osqp_solver), intent(inout) :: solver
      if (c_associated(solver%handle)) call c_of_api_destroy(solver%handle)
      solver%handle = c_null_ptr
      solver%initialized = .false.
      solver%model = osqp_model()
      solver%settings = osqp_settings()
   end subroutine osqp_cleanup

   subroutine osqp_solve_solver(solver, solution, status)
      type(osqp_solver), intent(inout) :: solver
      type(osqp_solution), intent(out) :: solution
      integer(osqp_int), intent(out), optional :: status
      real(c_double), allocatable, target :: x(:), y(:), pc(:), dc(:)
      type(of_info_c) :: ci
      character(kind=c_char), target :: buffer(128)
      integer(c_int) :: rc, rc2, nchar

      if (.not. solver%initialized .or. .not. c_associated(solver%handle)) then
         solution%call_status = osqp_workspace_not_init_error
         solution%status = "solver is not initialized"
         if (present(status)) status = solution%call_status
         return
      end if
      rc = c_of_api_solve(solver%handle)
      allocate(x(solver%model%n), y(solver%model%m), pc(solver%model%m), dc(solver%model%n))
      x = 0.0_dp; y = 0.0_dp; pc = 0.0_dp; dc = 0.0_dp
      rc2 = c_of_api_get_solution(solver%handle, pointer_or_null_real(x), pointer_or_null_real(y), &
         pointer_or_null_real(pc), pointer_or_null_real(dc), ci)
      if (rc2 == 0) then
         allocate(solution%x(size(x)), source=real(x,dp))
         allocate(solution%y(size(y)), source=real(y,dp))
         allocate(solution%prim_inf_cert(size(pc)), source=real(pc,dp))
         allocate(solution%dual_inf_cert(size(dc)), source=real(dc,dp))
         call info_from_c(ci, solution%info)
         nchar = c_of_api_status_message(solver%handle, c_loc(buffer(1)), int(size(buffer),c_int))
         solution%status = c_buffer_to_string(buffer, nchar)
      else
         solution%status = "solution could not be retrieved"
      end if
      solution%call_status = merge(rc, rc2, rc /= 0)
      if (present(status)) status = solution%call_status
   end subroutine osqp_solve_solver

   subroutine osqp_update(solver, status, q, l, u, px, px_idx, ax, ax_idx)
      type(osqp_solver), intent(inout) :: solver
      integer(osqp_int), intent(out) :: status
      real(dp), intent(in), optional :: q(:), l(:), u(:), px(:), ax(:)
      integer(osqp_int), intent(in), optional :: px_idx(:), ax_idx(:)
      real(c_double), allocatable, target :: cq(:), cl(:), cu(:), cpx(:), cax(:)
      integer(c_int), allocatable, target :: cpi(:), cai(:)
      type(c_ptr) :: qp, lp, up, pxp, pip, axp, aip
      integer(c_int) :: rc
      real(dp), allocatable :: test_l(:), test_u(:)

      status = osqp_invalid_argument
      if (.not. solver%initialized) return
      if (present(px_idx) .and. .not. present(px)) return
      if (present(ax_idx) .and. .not. present(ax)) return
      if (present(q)) then
         if (size(q) /= solver%model%n) return
         allocate(cq(size(q)), source=real(q,c_double))
      else
         allocate(cq(0))
      end if
      if (present(l)) then
         if (size(l) /= solver%model%m) return
         allocate(cl(size(l)), source=real(l,c_double))
      else
         allocate(cl(0))
      end if
      if (present(u)) then
         if (size(u) /= solver%model%m) return
         allocate(cu(size(u)), source=real(u,c_double))
      else
         allocate(cu(0))
      end if
      allocate(test_l(solver%model%m), source=solver%model%l)
      allocate(test_u(solver%model%m), source=solver%model%u)
      if (present(l)) test_l = l
      if (present(u)) test_u = u
      if (any(test_l > test_u)) return

      qp = pointer_or_null_real(cq)
      lp = pointer_or_null_real(cl)
      up = pointer_or_null_real(cu)
      rc = c_of_api_update_data_vec(solver%handle, qp, bool_int(present(q)), &
         lp, bool_int(present(l)), up, bool_int(present(u)))
      if (rc /= 0) then
         status = rc
         return
      end if
      if (present(q)) solver%model%q = q
      if (present(l)) solver%model%l = l
      if (present(u)) solver%model%u = u

      if (.not. present(px) .and. .not. present(ax)) then
         status = osqp_no_error
         return
      end if
      call prepare_matrix_update(px, px_idx, solver%model%p%value, cpx, cpi, status)
      if (status /= osqp_no_error) return
      call prepare_matrix_update(ax, ax_idx, solver%model%a%value, cax, cai, status)
      if (status /= osqp_no_error) return
      pxp = pointer_or_null_real(cpx)
      pip = pointer_or_null_int(cpi)
      axp = pointer_or_null_real(cax)
      aip = pointer_or_null_int(cai)
      rc = c_of_api_update_data_mat(solver%handle, pxp, pip, size(cpx), bool_int(present(px)), &
         axp, aip, size(cax), bool_int(present(ax)))
      status = rc
      if (rc /= 0) return
      if (present(px)) call apply_matrix_update(solver%model%p%value, px, px_idx)
      if (present(ax)) call apply_matrix_update(solver%model%a%value, ax, ax_idx)
   end subroutine osqp_update

   subroutine osqp_warm_start(solver, status, x, y)
      type(osqp_solver), intent(inout) :: solver
      integer(osqp_int), intent(out) :: status
      real(dp), intent(in), optional :: x(:), y(:)
      real(c_double), allocatable, target :: cx(:), cy(:)
      if (.not. solver%initialized) then
         status = osqp_workspace_not_init_error
         return
      end if
      if (.not. present(x) .and. .not. present(y)) then
         status = osqp_invalid_argument
         return
      end if
      if (present(x)) then
         if (size(x) /= solver%model%n) then
            status = osqp_invalid_argument
            return
         end if
         allocate(cx(size(x)), source=real(x,c_double))
      else
         allocate(cx(0))
      end if
      if (present(y)) then
         if (size(y) /= solver%model%m) then
            status = osqp_invalid_argument
            return
         end if
         allocate(cy(size(y)), source=real(y,c_double))
      else
         allocate(cy(0))
      end if
      status = c_of_api_warm_start(solver%handle, pointer_or_null_real(cx), bool_int(present(x)), &
         pointer_or_null_real(cy), bool_int(present(y)))
   end subroutine osqp_warm_start

   subroutine osqp_cold_start(solver, status)
      type(osqp_solver), intent(inout) :: solver
      integer(osqp_int), intent(out) :: status
      if (.not. solver%initialized) then
         status = osqp_workspace_not_init_error
      else
         status = c_of_api_cold_start(solver%handle)
      end if
   end subroutine osqp_cold_start

   subroutine osqp_update_settings(solver, settings, status)
      type(osqp_solver), intent(inout) :: solver
      type(osqp_settings), intent(in) :: settings
      integer(osqp_int), intent(out) :: status
      type(of_settings_c) :: cs
      if (.not. solver%initialized) then
         status = osqp_workspace_not_init_error
         return
      end if
      if (.not. settings%valid()) then
         status = osqp_invalid_argument
         return
      end if
      cs = settings_to_c(settings)
      status = c_of_api_update_settings(solver%handle, cs)
      if (status == 0) then
         status = c_of_api_get_settings(solver%handle, cs)
         if (status == 0) solver%settings = settings_from_c(cs)
      end if
   end subroutine osqp_update_settings

   subroutine osqp_get_settings(solver, settings, status)
      type(osqp_solver), intent(in) :: solver
      type(osqp_settings), intent(out) :: settings
      integer(osqp_int), intent(out), optional :: status
      type(of_settings_c) :: cs
      integer(c_int) :: rc
      if (.not. solver%initialized) then
         settings = osqp_settings()
         if (present(status)) status = osqp_workspace_not_init_error
         return
      end if
      rc = c_of_api_get_settings(solver%handle, cs)
      if (rc == 0) then
         settings = settings_from_c(cs)
      else
         settings = solver%settings
      end if
      if (present(status)) status = rc
   end subroutine osqp_get_settings

   subroutine osqp_get_dimensions(solver, n, m, status)
      type(osqp_solver), intent(in) :: solver
      integer(osqp_int), intent(out) :: n, m
      integer(osqp_int), intent(out), optional :: status
      integer(c_int) :: cn, cm, rc
      if (.not. solver%initialized) then
         n = 0; m = 0
         if (present(status)) status = osqp_workspace_not_init_error
         return
      end if
      rc = c_of_api_get_dimensions(solver%handle, cn, cm)
      n = cn; m = cm
      if (present(status)) status = rc
   end subroutine osqp_get_dimensions

   subroutine model_to_c_arrays(model, pp, pi, px, ap, ai, ax, q, l, u)
      type(osqp_model), intent(in) :: model
      integer(c_int), allocatable, target, intent(out) :: pp(:), pi(:), ap(:), ai(:)
      real(c_double), allocatable, target, intent(out) :: px(:), ax(:), q(:), l(:), u(:)
      allocate(pp(size(model%p%col_ptr)), source=int(model%p%col_ptr-1,c_int))
      allocate(pi(size(model%p%row_index)), source=int(model%p%row_index-1,c_int))
      allocate(px(size(model%p%value)), source=real(model%p%value,c_double))
      allocate(ap(size(model%a%col_ptr)), source=int(model%a%col_ptr-1,c_int))
      allocate(ai(size(model%a%row_index)), source=int(model%a%row_index-1,c_int))
      allocate(ax(size(model%a%value)), source=real(model%a%value,c_double))
      allocate(q(size(model%q)), source=real(model%q,c_double))
      allocate(l(size(model%l)), source=real(model%l,c_double))
      allocate(u(size(model%u)), source=real(model%u,c_double))
   end subroutine model_to_c_arrays

   subroutine prepare_matrix_update(values, indices, current, cvalues, cindices, status)
      real(dp), intent(in), optional :: values(:)
      integer(osqp_int), intent(in), optional :: indices(:)
      real(dp), intent(in) :: current(:)
      real(c_double), allocatable, target, intent(out) :: cvalues(:)
      integer(c_int), allocatable, target, intent(out) :: cindices(:)
      integer(osqp_int), intent(out) :: status
      status = osqp_no_error
      if (.not. present(values)) then
         allocate(cvalues(0), cindices(0))
         return
      end if
      if (present(indices)) then
         if (size(indices) /= size(values)) then
            status = osqp_invalid_argument
            allocate(cvalues(0), cindices(0))
            return
         end if
         if (size(indices) > 0) then
            if (any(indices < 1) .or. any(indices > size(current))) then
               status = osqp_invalid_argument
               allocate(cvalues(0), cindices(0))
               return
            end if
         end if
         allocate(cindices(size(indices)), source=int(indices-1,c_int))
      else
         if (size(values) /= size(current)) then
            status = osqp_invalid_argument
            allocate(cvalues(0), cindices(0))
            return
         end if
         allocate(cindices(0))
      end if
      allocate(cvalues(size(values)), source=real(values,c_double))
   end subroutine prepare_matrix_update

   subroutine apply_matrix_update(current, values, indices)
      real(dp), intent(inout) :: current(:)
      real(dp), intent(in) :: values(:)
      integer(osqp_int), intent(in), optional :: indices(:)
      integer :: k
      if (present(indices)) then
         do k = 1, size(values)
            current(indices(k)) = values(k)
         end do
      else
         current = values
      end if
   end subroutine apply_matrix_update

   pure type(of_settings_c) function settings_to_c(s) result(c)
      type(osqp_settings), intent(in) :: s
      c%device = s%device
      c%linsys_solver = s%linsys_solver
      c%allocate_solution = bool_int(s%allocate_solution)
      c%verbose = bool_int(s%verbose)
      c%profiler_level = s%profiler_level
      c%warm_starting = bool_int(s%warm_starting)
      c%scaling = s%scaling
      c%polishing = bool_int(s%polishing)
      c%rho = s%rho
      c%rho_is_vec = bool_int(s%rho_is_vec)
      c%sigma = s%sigma
      c%alpha = s%alpha
      c%cg_max_iter = s%cg_max_iter
      c%cg_tol_reduction = s%cg_tol_reduction
      c%cg_tol_fraction = s%cg_tol_fraction
      c%cg_precond = s%cg_precond
      c%adaptive_rho = s%adaptive_rho
      c%adaptive_rho_interval = s%adaptive_rho_interval
      c%adaptive_rho_fraction = s%adaptive_rho_fraction
      c%adaptive_rho_tolerance = s%adaptive_rho_tolerance
      c%max_iter = s%max_iter
      c%eps_abs = s%eps_abs
      c%eps_rel = s%eps_rel
      c%eps_prim_inf = s%eps_prim_inf
      c%eps_dual_inf = s%eps_dual_inf
      c%scaled_termination = bool_int(s%scaled_termination)
      c%check_termination = s%check_termination
      c%check_dualgap = bool_int(s%check_dualgap)
      c%time_limit = s%time_limit
      c%delta = s%delta
      c%polish_refine_iter = s%polish_refine_iter
   end function settings_to_c

   pure type(osqp_settings) function settings_from_c(c) result(s)
      type(of_settings_c), intent(in) :: c
      s%device = c%device
      s%linsys_solver = c%linsys_solver
      s%allocate_solution = c%allocate_solution /= 0
      s%verbose = c%verbose /= 0
      s%profiler_level = c%profiler_level
      s%warm_starting = c%warm_starting /= 0
      s%scaling = c%scaling
      s%polishing = c%polishing /= 0
      s%rho = c%rho
      s%rho_is_vec = c%rho_is_vec /= 0
      s%sigma = c%sigma
      s%alpha = c%alpha
      s%cg_max_iter = c%cg_max_iter
      s%cg_tol_reduction = c%cg_tol_reduction
      s%cg_tol_fraction = c%cg_tol_fraction
      s%cg_precond = c%cg_precond
      s%adaptive_rho = c%adaptive_rho
      s%adaptive_rho_interval = c%adaptive_rho_interval
      s%adaptive_rho_fraction = c%adaptive_rho_fraction
      s%adaptive_rho_tolerance = c%adaptive_rho_tolerance
      s%max_iter = c%max_iter
      s%eps_abs = c%eps_abs
      s%eps_rel = c%eps_rel
      s%eps_prim_inf = c%eps_prim_inf
      s%eps_dual_inf = c%eps_dual_inf
      s%scaled_termination = c%scaled_termination /= 0
      s%check_termination = c%check_termination
      s%check_dualgap = c%check_dualgap /= 0
      s%time_limit = c%time_limit
      s%delta = c%delta
      s%polish_refine_iter = c%polish_refine_iter
   end function settings_from_c

   pure subroutine info_from_c(c, info)
      type(of_info_c), intent(in) :: c
      type(osqp_info), intent(out) :: info
      info%status_val = c%status_val
      info%status_polish = c%status_polish
      info%obj_val = c%obj_val
      info%dual_obj_val = c%dual_obj_val
      info%prim_res = c%prim_res
      info%dual_res = c%dual_res
      info%duality_gap = c%duality_gap
      info%iter = c%iter
      info%rho_updates = c%rho_updates
      info%rho_estimate = c%rho_estimate
      info%setup_time = c%setup_time
      info%solve_time = c%solve_time
      info%update_time = c%update_time
      info%polish_time = c%polish_time
      info%run_time = c%run_time
      info%primdual_int = c%primdual_int
      info%rel_kkt_error = c%rel_kkt_error
   end subroutine info_from_c

   pure integer(c_int) function bool_int(x) result(i)
      logical, intent(in) :: x
      i = merge(1_c_int, 0_c_int, x)
   end function bool_int

   function pointer_or_null_real(x) result(p)
      real(c_double), allocatable, target, intent(inout) :: x(:)
      type(c_ptr) :: p
      if (size(x) > 0) then
         p = c_loc(x(1))
      else
         p = c_null_ptr
      end if
   end function pointer_or_null_real

   function pointer_or_null_int(x) result(p)
      integer(c_int), allocatable, target, intent(inout) :: x(:)
      type(c_ptr) :: p
      if (size(x) > 0) then
         p = c_loc(x(1))
      else
         p = c_null_ptr
      end if
   end function pointer_or_null_int

   function c_buffer_to_string(buffer, reported_length) result(text)
      character(kind=c_char), intent(in) :: buffer(:)
      integer(c_int), intent(in) :: reported_length
      character(len=:), allocatable :: text
      integer :: i, n
      n = min(int(reported_length), size(buffer))
      do i = 1, n
         if (buffer(i) == c_null_char) then
            n = i - 1
            exit
         end if
      end do
      allocate(character(len=max(0,n)) :: text)
      do i = 1, n
         text(i:i) = achar(iachar(buffer(i)))
      end do
   end function c_buffer_to_string

end module osqp_solver_api
