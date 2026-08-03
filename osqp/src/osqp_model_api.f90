! SPDX-License-Identifier: Apache-2.0
module osqp_model_api
   use osqp_kinds, only : dp, osqp_int
   use osqp_constants
   use osqp_sparse
   use osqp_types
   use osqp_solver_api
   implicit none
   private

   public :: osqp_model_from_dense, osqp_model_from_sparse, solve_osqp

   interface solve_osqp
      module procedure solve_osqp_model
      module procedure solve_osqp_arrays
   end interface solve_osqp

contains

   subroutine osqp_model_from_dense(model, q, status, p, a, l, u, tolerance)
      type(osqp_model), intent(out) :: model
      real(dp), intent(in) :: q(:)
      integer(osqp_int), intent(out) :: status
      real(dp), intent(in), optional :: p(:,:), a(:,:), l(:), u(:)
      real(dp), intent(in), optional :: tolerance
      integer :: n, m
      real(dp) :: tol

      status = osqp_invalid_argument
      n = size(q)
      if (n < 1) return
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      model%n = n
      allocate(model%q(n), source=q)
      if (present(p)) then
         if (size(p,1) /= n .or. size(p,2) /= n) return
         model%p = osqp_csc_from_dense(p, upper_only=.true., tolerance=tol)
      else
         model%p = osqp_empty_matrix(n, n)
      end if
      if (present(a)) then
         if (size(a,2) /= n) return
         m = size(a,1)
         model%m = m
         model%a = osqp_csc_from_dense(a, tolerance=tol)
         allocate(model%l(m), model%u(m))
         model%l = -osqp_infinity_default
         model%u = osqp_infinity_default
         if (present(l)) then
            if (size(l) /= m) return
            model%l = l
         end if
         if (present(u)) then
            if (size(u) /= m) return
            model%u = u
         end if
      else
         if (present(l) .or. present(u)) return
         model%m = 0
         model%a = osqp_empty_matrix(0, n)
         allocate(model%l(0), model%u(0))
      end if
      if (.not. model%valid()) return
      status = osqp_no_error
   end subroutine osqp_model_from_dense

   subroutine osqp_model_from_sparse(model, p, q, a, l, u, status)
      type(osqp_model), intent(out) :: model
      type(osqp_sparse_matrix), intent(in) :: p, a
      real(dp), intent(in) :: q(:), l(:), u(:)
      integer(osqp_int), intent(out) :: status
      status = osqp_invalid_argument
      model%n = size(q)
      model%m = size(l)
      if (size(u) /= model%m) return
      model%p = p
      model%a = a
      allocate(model%q(size(q)), source=q)
      allocate(model%l(size(l)), source=l)
      allocate(model%u(size(u)), source=u)
      if (.not. model%valid()) return
      status = osqp_no_error
   end subroutine osqp_model_from_sparse

   subroutine solve_osqp_model(model, solution, status, settings, warm_x, warm_y)
      type(osqp_model), intent(in) :: model
      type(osqp_solution), intent(out) :: solution
      integer(osqp_int), intent(out), optional :: status
      type(osqp_settings), intent(in), optional :: settings
      real(dp), intent(in), optional :: warm_x(:), warm_y(:)
      type(osqp_solver) :: solver
      integer(osqp_int) :: rc

      call osqp_setup(solver, model, rc, settings)
      if (rc /= osqp_no_error) then
         solution%call_status = rc
         solution%status = osqp_error_string(rc)
         if (present(status)) status = rc
         return
      end if
      if (present(warm_x) .or. present(warm_y)) then
         call osqp_warm_start(solver, rc, warm_x, warm_y)
         if (rc /= osqp_no_error) then
            solution%call_status = rc
            solution%status = osqp_error_string(rc)
            call osqp_cleanup(solver)
            if (present(status)) status = rc
            return
         end if
      end if
      call osqp_solve_solver(solver, solution, rc)
      call osqp_cleanup(solver)
      if (present(status)) status = rc
   end subroutine solve_osqp_model

   subroutine solve_osqp_arrays(q, solution, status, p, a, l, u, settings, tolerance, warm_x, warm_y)
      real(dp), intent(in) :: q(:)
      type(osqp_solution), intent(out) :: solution
      integer(osqp_int), intent(out), optional :: status
      real(dp), intent(in), optional :: p(:,:), a(:,:), l(:), u(:), tolerance
      type(osqp_settings), intent(in), optional :: settings
      real(dp), intent(in), optional :: warm_x(:), warm_y(:)
      type(osqp_model) :: model
      integer(osqp_int) :: rc

      call osqp_model_from_dense(model, q, rc, p, a, l, u, tolerance)
      if (rc /= osqp_no_error) then
         solution%call_status = rc
         solution%status = "invalid model dimensions or data"
         if (present(status)) status = rc
         return
      end if
      call solve_osqp_model(model, solution, rc, settings, warm_x, warm_y)
      if (present(status)) status = rc
   end subroutine solve_osqp_arrays

end module osqp_model_api
