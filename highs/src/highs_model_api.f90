! SPDX-License-Identifier: GPL-2.0-or-later
module highs_model_api
   use highs_kinds, only : dp, highs_int
   use highs_constants
   use highs_sparse
   use highs_types
   use highs_solver_api
   implicit none
   private

   public :: highs_model_from_dense, highs_solve

   interface highs_solve
      module procedure highs_solve_model
      module procedure highs_solve_arrays
   end interface highs_solve

contains

   subroutine highs_model_from_dense(model, linear, lower, upper, status, a, lhs, rhs, &
         vartype, maximum, offset, q)
      type(highs_model), intent(out) :: model
      real(dp), intent(in) :: linear(:), lower(:), upper(:)
      integer(highs_int), intent(out) :: status
      real(dp), intent(in), optional :: a(:,:), lhs(:), rhs(:), q(:,:)
      integer(highs_int), intent(in), optional :: vartype(:)
      logical, intent(in), optional :: maximum
      real(dp), intent(in), optional :: offset
      integer :: n, m

      status = highs_invalid_argument
      n = size(linear)
      if (size(lower) /= n .or. size(upper) /= n) return
      model%num_col = n
      allocate(model%col_cost(n), source=linear)
      allocate(model%col_lower(n), source=lower)
      allocate(model%col_upper(n), source=upper)
      allocate(model%integrality(n), source=highs_var_continuous)
      if (present(vartype)) then
         if (size(vartype) /= n .or. any(vartype < 0) .or. any(vartype > 4)) return
         model%integrality = vartype
      end if
      if (present(maximum)) then
         model%sense = merge(highs_maximize, highs_minimize, maximum)
      else
         model%sense = highs_minimize
      end if
      model%offset = 0.0_dp
      if (present(offset)) model%offset = offset

      if (present(a)) then
         m = size(a,1)
         if (size(a,2) /= n) return
         if (.not. present(lhs) .or. .not. present(rhs)) return
         if (size(lhs) /= m .or. size(rhs) /= m) return
         model%num_row = m
         allocate(model%row_lower(m), source=lhs)
         allocate(model%row_upper(m), source=rhs)
         model%a = highs_csc_from_dense(a)
      else
         if (present(lhs) .or. present(rhs)) return
         model%num_row = 0
         allocate(model%row_lower(0), model%row_upper(0))
         model%a = highs_empty_matrix(0, n, highs_matrix_colwise)
      end if

      if (present(q)) then
         if (size(q,1) /= n .or. size(q,2) /= n) return
         if (maxval(abs(q - transpose(q))) > 100.0_dp * epsilon(1.0_dp) * &
             max(1.0_dp, maxval(abs(q)))) return
         model%q = highs_hessian_from_dense(q)
         model%has_hessian = .true.
      else
         model%q = highs_empty_matrix(n, n, highs_hessian_triangular)
         model%has_hessian = .false.
      end if
      if (.not. model%valid()) return
      status = highs_status_ok
   end subroutine highs_model_from_dense

   subroutine highs_solve_model(model, solution, status, control, start)
      type(highs_model), intent(in) :: model
      type(highs_solution), intent(out) :: solution
      integer(highs_int), intent(out), optional :: status
      type(highs_control), intent(in), optional :: control
      real(dp), intent(in), optional :: start(:)
      type(highs_solver) :: solver
      type(highs_control) :: ctrl
      integer(highs_int) :: rc

      call highs_new_solver(solver, rc)
      if (rc /= highs_status_ok) then
         solution%call_status = rc
         solution%status_message = "HiGHS backend unavailable: " // highs_backend_error()
         if (present(status)) status = rc
         return
      end if
      call highs_pass_model(solver, model, rc)
      if (rc == highs_status_error .or. rc == highs_invalid_argument) then
         solution%call_status = rc
         solution%status_message = "model could not be passed to HiGHS"
         call highs_destroy_solver(solver)
         if (present(status)) status = rc
         return
      end if
      ctrl = highs_control()
      if (present(control)) ctrl = control
      call highs_set_control(solver, ctrl, rc)
      if (present(start)) then
         call highs_set_start(solver, start, rc)
      end if
      call highs_run(solver, rc)
      call highs_get_solution(solver, solution)
      solution%call_status = rc
      call highs_destroy_solver(solver)
      if (present(status)) status = rc
   end subroutine highs_solve_model

   subroutine highs_solve_arrays(linear, lower, upper, solution, status, a, lhs, rhs, &
         vartype, maximum, offset, q, control, start)
      real(dp), intent(in) :: linear(:), lower(:), upper(:)
      type(highs_solution), intent(out) :: solution
      integer(highs_int), intent(out), optional :: status
      real(dp), intent(in), optional :: a(:,:), lhs(:), rhs(:), q(:,:), start(:)
      integer(highs_int), intent(in), optional :: vartype(:)
      logical, intent(in), optional :: maximum
      real(dp), intent(in), optional :: offset
      type(highs_control), intent(in), optional :: control
      type(highs_model) :: model
      integer(highs_int) :: rc

      call highs_model_from_dense(model, linear, lower, upper, rc, a, lhs, rhs, &
         vartype, maximum, offset, q)
      if (rc /= highs_status_ok) then
         solution%call_status = rc
         solution%status_message = "invalid model dimensions or data"
         if (present(status)) status = rc
         return
      end if
      call highs_solve_model(model, solution, rc, control, start)
      if (present(status)) status = rc
   end subroutine highs_solve_arrays

end module highs_model_api
