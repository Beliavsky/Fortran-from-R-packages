! SPDX-License-Identifier: Apache-2.0
module clarabel_solver
   use, intrinsic :: iso_c_binding
   use clarabel_kinds, only : dp
   use clarabel_sparse, only : csc_matrix
   use clarabel_types
   use clarabel_c_api
   implicit none
   private

   integer, parameter :: error_capacity = 1024
   integer, parameter, public :: clarabel_backend_unavailable = -100

   type, public :: clarabel_solver_type
      private
      type(c_ptr) :: handle = c_null_ptr
      integer(c_size_t) :: n = 0_c_size_t
      integer(c_size_t) :: m = 0_c_size_t
      integer(c_size_t) :: p_nnz = 0_c_size_t
      integer(c_size_t) :: a_nnz = 0_c_size_t
      integer(c_size_t), allocatable :: p_colptr(:), p_rowind(:)
      integer(c_size_t), allocatable :: a_colptr(:), a_rowind(:)
   contains
      procedure, public :: initialize => solver_initialize
      procedure, public :: solve => solver_solve
      procedure, public :: update => solver_update
      procedure, public :: is_update_allowed => solver_is_update_allowed
      procedure, public :: release => solver_release
      procedure, public :: valid => solver_valid
      final :: solver_finalize
   end type clarabel_solver_type

   public :: clarabel_solve_problem

contains

   subroutine clarabel_solve_problem(p, q, a, b, cones, solution, settings, code, message)
      type(csc_matrix), intent(in) :: p, a
      real(dp), intent(in) :: q(:), b(:)
      type(clarabel_cone), intent(in) :: cones(:)
      type(clarabel_solution), intent(out) :: solution
      type(clarabel_settings), intent(in), optional :: settings
      integer, intent(out), optional :: code
      character(len=:), allocatable, intent(out), optional :: message
      type(clarabel_solver_type) :: solver
      integer :: ierr
      character(len=:), allocatable :: text

      call solver%initialize(p, q, a, b, cones, settings, ierr, text)
      if (ierr == 0) call solver%solve(solution, ierr, text)
      if (present(code)) code = ierr
      if (present(message)) message = text
   end subroutine clarabel_solve_problem

   subroutine solver_initialize(self, p, q, a, b, cones, settings, code, message)
      class(clarabel_solver_type), intent(inout) :: self
      type(csc_matrix), intent(in) :: p, a
      real(dp), intent(in) :: q(:), b(:)
      type(clarabel_cone), intent(in) :: cones(:)
      type(clarabel_settings), intent(in), optional :: settings
      integer, intent(out), optional :: code
      character(len=:), allocatable, intent(out), optional :: message

      type(clarabel_csc_c) :: pc, ac
      type(clarabel_settings_c) :: sc
      type(clarabel_cone_c), allocatable, target :: cc(:)
      integer(c_size_t), allocatable, target :: p_col(:), p_row(:), a_col(:), a_row(:)
      real(dp), allocatable, target :: p_val(:), a_val(:), q_work(:), b_work(:), alpha(:)
      character(c_char) :: error(error_capacity)
      type(c_ptr) :: q_ptr, b_ptr, cone_ptr
      integer(c_int32_t) :: rc
      integer :: ierr
      character(len=:), allocatable :: text

      call self%release()
      call validate_problem(p, q, a, b, cones, ierr, text)
      if (ierr /= 0) then
         call set_outputs(code, message, ierr, text)
         return
      end if

      p_col = p%colptr
      p_row = p%rowind
      p_val = p%values
      a_col = a%colptr
      a_row = a%rowind
      a_val = a%values
      q_work = q
      b_work = b

      call make_csc_descriptor(p, p_col, p_row, p_val, pc)
      call make_csc_descriptor(a, a_col, a_row, a_val, ac)
      call pack_cones(cones, cc, alpha)
      if (size(q_work) > 0) then
         q_ptr = c_loc(q_work(1))
      else
         q_ptr = c_null_ptr
      end if
      if (size(b_work) > 0) then
         b_ptr = c_loc(b_work(1))
      else
         b_ptr = c_null_ptr
      end if
      if (size(cc) > 0) then
         cone_ptr = c_loc(cc(1))
      else
         cone_ptr = c_null_ptr
      end if

      if (present(settings)) then
         call settings_to_c(settings, sc)
      else
         call settings_to_c(default_clarabel_settings(), sc)
      end if
      error = c_null_char
      rc = c_clarabel_solver_create(pc, q_ptr, int(size(q_work), c_size_t), ac, b_ptr, &
                                    int(size(b_work), c_size_t), cone_ptr, int(size(cc), c_size_t), &
                                    sc, self%handle, error, int(error_capacity, c_size_t))
      if (rc /= 0_c_int32_t .or. .not. c_associated(self%handle)) then
         self%handle = c_null_ptr
         ierr = int(rc)
         if (ierr == 0) ierr = -1
         text = c_buffer_to_string(error)
         if (len(text) == 0) text = "Clarabel backend could not create the solver"
         call set_outputs(code, message, ierr, text)
         return
      end if

      self%n = p%ncols
      self%m = a%nrows
      self%p_nnz = int(size(p%values), c_size_t)
      self%a_nnz = int(size(a%values), c_size_t)
      self%p_colptr = p%colptr
      self%p_rowind = p%rowind
      self%a_colptr = a%colptr
      self%a_rowind = a%rowind
      call set_outputs(code, message, 0, "")
   end subroutine solver_initialize

   subroutine solver_solve(self, solution, code, message)
      class(clarabel_solver_type), intent(inout) :: self
      type(clarabel_solution), intent(out) :: solution
      integer, intent(out), optional :: code
      character(len=:), allocatable, intent(out), optional :: message
      real(dp), allocatable, target :: x(:), z(:), s(:)
      type(clarabel_result_c) :: result
      character(c_char) :: error(error_capacity)
      type(c_ptr) :: x_ptr, z_ptr, s_ptr
      integer(c_int32_t) :: rc
      integer :: ierr
      character(len=:), allocatable :: text

      if (.not. c_associated(self%handle)) then
         call set_outputs(code, message, -1, "solver is not initialized")
         return
      end if
      allocate(x(int(self%n)), z(int(self%m)), s(int(self%m)))
      x = 0.0_dp
      z = 0.0_dp
      s = 0.0_dp
      x_ptr = pointer_or_null(x)
      z_ptr = pointer_or_null(z)
      s_ptr = pointer_or_null(s)
      error = c_null_char
      rc = c_clarabel_solver_solve(self%handle, x_ptr, self%n, z_ptr, self%m, s_ptr, self%m, &
                                   result, error, int(error_capacity, c_size_t))
      if (rc /= 0_c_int32_t) then
         ierr = int(rc)
         text = c_buffer_to_string(error)
         if (len(text) == 0) text = "Clarabel backend solve failed"
         call set_outputs(code, message, ierr, text)
         return
      end if
      call move_alloc(x, solution%x)
      call move_alloc(z, solution%z)
      call move_alloc(s, solution%s)
      call result_from_c(result, solution)
      call set_outputs(code, message, 0, "")
   end subroutine solver_solve

   subroutine solver_update(self, p, a, q, b, code, message)
      class(clarabel_solver_type), intent(inout) :: self
      type(csc_matrix), intent(in), optional :: p, a
      real(dp), intent(in), optional :: q(:), b(:)
      integer, intent(out), optional :: code
      character(len=:), allocatable, intent(out), optional :: message
      real(dp), allocatable, target :: pv(:), av(:), qv(:), bv(:)
      type(c_ptr) :: pp, ap, qp, bp
      integer(c_size_t) :: pn, an, qn, bn
      character(c_char) :: error(error_capacity)
      integer(c_int32_t) :: rc
      integer :: ierr
      character(len=:), allocatable :: text

      if (.not. c_associated(self%handle)) then
         call set_outputs(code, message, -1, "solver is not initialized")
         return
      end if
      pp = c_null_ptr; ap = c_null_ptr; qp = c_null_ptr; bp = c_null_ptr
      pn = 0_c_size_t; an = 0_c_size_t; qn = 0_c_size_t; bn = 0_c_size_t
      if (present(p)) then
         if (p%nrows /= self%n .or. p%ncols /= self%n .or. &
             int(size(p%values), c_size_t) /= self%p_nnz) then
            call set_outputs(code, message, -2, "updated P has incompatible dimensions or nonzero count")
            return
         end if
         if (any(p%colptr /= self%p_colptr) .or. any(p%rowind /= self%p_rowind)) then
            call set_outputs(code, message, -2, "updated P must preserve the original sparsity pattern")
            return
         end if
         pv = p%values
         pn = int(size(pv), c_size_t)
         pp = pointer_or_null(pv)
      end if
      if (present(a)) then
         if (a%nrows /= self%m .or. a%ncols /= self%n .or. &
             int(size(a%values), c_size_t) /= self%a_nnz) then
            call set_outputs(code, message, -3, "updated A has incompatible dimensions or nonzero count")
            return
         end if
         if (any(a%colptr /= self%a_colptr) .or. any(a%rowind /= self%a_rowind)) then
            call set_outputs(code, message, -3, "updated A must preserve the original sparsity pattern")
            return
         end if
         av = a%values
         an = int(size(av), c_size_t)
         ap = pointer_or_null(av)
      end if
      if (present(q)) then
         if (int(size(q), c_size_t) /= self%n) then
            call set_outputs(code, message, -4, "updated q has the wrong length")
            return
         end if
         qv = q
         qn = int(size(qv), c_size_t)
         qp = pointer_or_null(qv)
      end if
      if (present(b)) then
         if (int(size(b), c_size_t) /= self%m) then
            call set_outputs(code, message, -5, "updated b has the wrong length")
            return
         end if
         bv = b
         bn = int(size(bv), c_size_t)
         bp = pointer_or_null(bv)
      end if
      error = c_null_char
      rc = c_clarabel_solver_update(self%handle, pp, pn, ap, an, qp, qn, bp, bn, error, &
                                    int(error_capacity, c_size_t))
      if (rc /= 0_c_int32_t) then
         ierr = int(rc)
         text = c_buffer_to_string(error)
         if (len(text) == 0) text = "Clarabel backend update failed"
         call set_outputs(code, message, ierr, text)
         return
      end if
      call set_outputs(code, message, 0, "")
   end subroutine solver_update

   logical function solver_is_update_allowed(self) result(allowed)
      class(clarabel_solver_type), intent(in) :: self
      if (c_associated(self%handle)) then
         allowed = c_clarabel_solver_is_update_allowed(self%handle) /= 0_c_int32_t
      else
         allowed = .false.
      end if
   end function solver_is_update_allowed

   logical function solver_valid(self) result(ok)
      class(clarabel_solver_type), intent(in) :: self
      ok = c_associated(self%handle)
   end function solver_valid

   subroutine solver_release(self)
      class(clarabel_solver_type), intent(inout) :: self
      if (c_associated(self%handle)) call c_clarabel_solver_free(self%handle)
      self%handle = c_null_ptr
      self%n = 0_c_size_t
      self%m = 0_c_size_t
      self%p_nnz = 0_c_size_t
      self%a_nnz = 0_c_size_t
      if (allocated(self%p_colptr)) deallocate(self%p_colptr)
      if (allocated(self%p_rowind)) deallocate(self%p_rowind)
      if (allocated(self%a_colptr)) deallocate(self%a_colptr)
      if (allocated(self%a_rowind)) deallocate(self%a_rowind)
   end subroutine solver_release

   subroutine solver_finalize(self)
      type(clarabel_solver_type), intent(inout) :: self
      call self%release()
   end subroutine solver_finalize

   subroutine validate_problem(p, q, a, b, cones, code, message)
      type(csc_matrix), intent(in) :: p, a
      real(dp), intent(in) :: q(:), b(:)
      type(clarabel_cone), intent(in) :: cones(:)
      integer, intent(out) :: code
      character(len=:), allocatable, intent(out) :: message
      logical :: ok
      integer :: i, j, k

      call p%validate(ok, message)
      if (.not. ok) then
         code = -10
         message = "P: " // message
         return
      end if
      call a%validate(ok, message)
      if (.not. ok) then
         code = -11
         message = "A: " // message
         return
      end if
      if (p%nrows /= p%ncols) then
         code = -12; message = "P must be square"; return
      end if
      do j = 1, int(p%ncols)
         do k = int(p%colptr(j)) + 1, int(p%colptr(j + 1))
            if (int(p%rowind(k)) + 1 > j) then
               code = -18; message = "P must contain only its upper triangle"; return
            end if
         end do
      end do
      if (a%ncols /= p%ncols) then
         code = -13; message = "A and P have incompatible column counts"; return
      end if
      if (size(q) /= int(p%ncols)) then
         code = -14; message = "q has the wrong length"; return
      end if
      if (size(b) /= int(a%nrows)) then
         code = -15; message = "b has the wrong length"; return
      end if
      if (cones_total_dimension(cones) /= a%nrows) then
         code = -16; message = "cone dimensions do not equal the row count of A"; return
      end if
      do i = 1, size(cones)
         call cones(i)%validate(ok, message)
         if (.not. ok) then
            code = -17
            message = "cone " // integer_string(i) // ": " // message
            return
         end if
      end do
      code = 0
      message = ""
   end subroutine validate_problem

   subroutine make_csc_descriptor(source, colptr, rowind, values, out)
      type(csc_matrix), intent(in) :: source
      integer(c_size_t), target, intent(in) :: colptr(:), rowind(:)
      real(dp), target, intent(in) :: values(:)
      type(clarabel_csc_c), intent(out) :: out
      out%nrows = source%nrows
      out%ncols = source%ncols
      out%nnz = int(size(values), c_size_t)
      out%colptr = c_loc(colptr(1))
      if (size(rowind) > 0) then
         out%rowind = c_loc(rowind(1))
         out%values = c_loc(values(1))
      else
         out%rowind = c_null_ptr
         out%values = c_null_ptr
      end if
   end subroutine make_csc_descriptor

   subroutine pack_cones(cones, packed, alpha)
      type(clarabel_cone), intent(in) :: cones(:)
      type(clarabel_cone_c), allocatable, target, intent(out) :: packed(:)
      real(dp), allocatable, target, intent(out) :: alpha(:)
      integer :: i, k, nalpha

      nalpha = 0
      do i = 1, size(cones)
         if (allocated(cones(i)%alpha)) nalpha = nalpha + size(cones(i)%alpha)
      end do
      allocate(packed(size(cones)), alpha(nalpha))
      k = 1
      do i = 1, size(cones)
         packed(i)%tag = cones(i)%kind
         packed(i)%dim = cones(i)%dim
         packed(i)%parameter = cones(i)%parameter
         if (allocated(cones(i)%alpha)) then
            packed(i)%alpha_len = int(size(cones(i)%alpha), c_size_t)
            alpha(k:k + size(cones(i)%alpha) - 1) = cones(i)%alpha
            packed(i)%alpha = c_loc(alpha(k))
            k = k + size(cones(i)%alpha)
         else
            packed(i)%alpha_len = 0_c_size_t
            packed(i)%alpha = c_null_ptr
         end if
      end do
   end subroutine pack_cones

   subroutine settings_to_c(source, target)
      type(clarabel_settings), intent(in) :: source
      type(clarabel_settings_c), intent(out) :: target
      target%max_iter = source%max_iter
      target%time_limit = source%time_limit
      target%verbose = logical_to_c(source%verbose)
      target%max_step_fraction = source%max_step_fraction
      target%tol_gap_abs = source%tol_gap_abs
      target%tol_gap_rel = source%tol_gap_rel
      target%tol_feas = source%tol_feas
      target%tol_infeas_abs = source%tol_infeas_abs
      target%tol_infeas_rel = source%tol_infeas_rel
      target%tol_ktratio = source%tol_ktratio
      target%reduced_tol_gap_abs = source%reduced_tol_gap_abs
      target%reduced_tol_gap_rel = source%reduced_tol_gap_rel
      target%reduced_tol_feas = source%reduced_tol_feas
      target%reduced_tol_infeas_abs = source%reduced_tol_infeas_abs
      target%reduced_tol_infeas_rel = source%reduced_tol_infeas_rel
      target%reduced_tol_ktratio = source%reduced_tol_ktratio
      target%equilibrate_enable = logical_to_c(source%equilibrate_enable)
      target%equilibrate_max_iter = source%equilibrate_max_iter
      target%equilibrate_min_scaling = source%equilibrate_min_scaling
      target%equilibrate_max_scaling = source%equilibrate_max_scaling
      target%linesearch_backtrack_step = source%linesearch_backtrack_step
      target%min_switch_step_length = source%min_switch_step_length
      target%min_terminate_step_length = source%min_terminate_step_length
      target%max_threads = source%max_threads
      target%direct_kkt_solver = logical_to_c(source%direct_kkt_solver)
      target%direct_solve_method = source%direct_solve_method
      target%static_regularization_enable = logical_to_c(source%static_regularization_enable)
      target%static_regularization_constant = source%static_regularization_constant
      target%static_regularization_proportional = source%static_regularization_proportional
      target%dynamic_regularization_enable = logical_to_c(source%dynamic_regularization_enable)
      target%dynamic_regularization_eps = source%dynamic_regularization_eps
      target%dynamic_regularization_delta = source%dynamic_regularization_delta
      target%iterative_refinement_enable = logical_to_c(source%iterative_refinement_enable)
      target%iterative_refinement_reltol = source%iterative_refinement_reltol
      target%iterative_refinement_abstol = source%iterative_refinement_abstol
      target%iterative_refinement_max_iter = source%iterative_refinement_max_iter
      target%iterative_refinement_stop_ratio = source%iterative_refinement_stop_ratio
      target%presolve_enable = logical_to_c(source%presolve_enable)
      target%input_sparse_dropzeros = logical_to_c(source%input_sparse_dropzeros)
      target%chordal_decomposition_enable = logical_to_c(source%chordal_decomposition_enable)
      target%chordal_decomposition_merge_method = source%chordal_decomposition_merge_method
      target%chordal_decomposition_compact = logical_to_c(source%chordal_decomposition_compact)
      target%chordal_decomposition_complete_dual = logical_to_c(source%chordal_decomposition_complete_dual)
   end subroutine settings_to_c

   subroutine result_from_c(source, target)
      type(clarabel_result_c), intent(in) :: source
      type(clarabel_solution), intent(inout) :: target
      target%status = source%status
      target%iterations = source%iterations
      target%obj_val = source%obj_val
      target%obj_val_dual = source%obj_val_dual
      target%solve_time = source%solve_time
      target%r_prim = source%r_prim
      target%r_dual = source%r_dual
      target%info%status = source%status
      target%info%iterations = source%iterations
      target%info%mu = source%mu
      target%info%sigma = source%sigma
      target%info%step_length = source%step_length
      target%info%cost_primal = source%cost_primal
      target%info%cost_dual = source%cost_dual
      target%info%res_primal = source%res_primal
      target%info%res_dual = source%res_dual
      target%info%res_primal_inf = source%res_primal_inf
      target%info%res_dual_inf = source%res_dual_inf
      target%info%gap_abs = source%gap_abs
      target%info%gap_rel = source%gap_rel
      target%info%ktratio = source%ktratio
      target%info%solve_time = source%solve_time
      target%info%linear_solver_threads = source%linear_solver_threads
      target%info%linear_solver_nnz_a = source%linear_solver_nnz_a
      target%info%linear_solver_nnz_l = source%linear_solver_nnz_l
   end subroutine result_from_c

   pure integer(c_int32_t) function logical_to_c(x) result(value)
      logical, intent(in) :: x
      if (x) then
         value = 1_c_int32_t
      else
         value = 0_c_int32_t
      end if
   end function logical_to_c

   function pointer_or_null(x) result(ptr)
      real(dp), target, intent(inout) :: x(:)
      type(c_ptr) :: ptr
      if (size(x) > 0) then
         ptr = c_loc(x(1))
      else
         ptr = c_null_ptr
      end if
   end function pointer_or_null

   function c_buffer_to_string(buffer) result(text)
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
   end function c_buffer_to_string

   subroutine set_outputs(code, message, value, text)
      integer, intent(out), optional :: code
      character(len=:), allocatable, intent(out), optional :: message
      integer, intent(in) :: value
      character(len=*), intent(in) :: text
      if (present(code)) code = value
      if (present(message)) message = text
   end subroutine set_outputs

   pure function integer_string(i) result(text)
      integer, intent(in) :: i
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      write(buffer, '(i0)') i
      text = trim(buffer)
   end function integer_string

end module clarabel_solver
