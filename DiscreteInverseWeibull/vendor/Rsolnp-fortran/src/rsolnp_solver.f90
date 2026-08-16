! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rsolnp_kinds, only : dp
   use rsolnp_types, only : solnp_problem, solnp_control, solnp_result, kkt_diagnostics, &
      solnp_success, solnp_max_iterations, solnp_invalid_problem, solnp_numerical_failure, &
      solnp_infeasible
   use rsolnp_problem, only : prepare_problem
   use rsolnp_evaluate, only : evaluate_objective, evaluate_gradient, evaluate_equalities, &
      evaluate_inequalities, evaluate_eq_jacobian, evaluate_ineq_jacobian, constraint_values, &
      constraint_jacobian
   use rsolnp_linalg, only : vnorm2, norm_inf, eye, outer_product, invert_matrix, &
      symmetrize, projected_gradient
   implicit none
   private

   public :: solnp, csolnp, kkt_diagnose

contains

   subroutine solnp(problem, result, control)
      type(solnp_problem), intent(in) :: problem
      type(solnp_result), intent(out) :: result
      type(solnp_control), intent(in), optional :: control
      call csolnp(problem, result, control)
   end subroutine solnp

   subroutine csolnp(problem, result, control)
      type(solnp_problem), intent(in) :: problem
      type(solnp_result), intent(out) :: result
      type(solnp_control), intent(in), optional :: control

      type(solnp_problem) :: prob
      type(solnp_control) :: ctrl
      real(dp), allocatable :: y(:), lower(:), upper(:), lambda(:), c(:), h_inv(:, :)
      real(dp), allocatable :: obj_hist(:), con_hist(:), step_hist(:), full_hess(:, :)
      real(dp) :: rho, f, previous_f, cnorm, previous_cnorm, rel_obj, step_norm
      real(dp) :: start_time, end_time
      integer :: n, m, nc, outer, inner_used, total_inner, stat
      logical :: ok, inv_ok
      character(len=160) :: msg

      call cpu_time(start_time)
      ctrl = solnp_control()
      if (present(control)) ctrl = control
      prob = problem
      call prepare_problem(prob, stat, msg)
      if (stat /= solnp_success) then
         result%convergence = solnp_invalid_problem
         result%message = trim(msg)
         call cpu_time(end_time)
         result%elapsed = end_time - start_time
         return
      end if
      if (ctrl%max_iter < 1 .or. ctrl%min_iter < 1 .or. ctrl%tol <= 0.0_dp .or. &
          ctrl%delta <= 0.0_dp) then
         result%convergence = solnp_invalid_problem
         result%message = 'invalid solver control values'
         call cpu_time(end_time)
         result%elapsed = end_time - start_time
         return
      end if

      n = prob%n
      m = prob%n_ineq
      nc = prob%n_eq + prob%n_ineq
      allocate(y(n + m), lower(n + m), upper(n + m), lambda(nc), c(nc))
      allocate(h_inv(n + m, n + m))
      allocate(obj_hist(ctrl%max_iter + 1), con_hist(ctrl%max_iter + 1), &
               step_hist(ctrl%max_iter + 1))
      y(m + 1:) = prob%start
      lower(m + 1:) = prob%lower
      upper(m + 1:) = prob%upper
      if (m > 0) then
         lower(1:m) = prob%ineq_lower
         upper(1:m) = prob%ineq_upper
         call initialize_slack(prob, prob%start, y(1:m), ok)
         if (.not. ok) then
            result%convergence = solnp_numerical_failure
            result%message = 'nonfinite initial inequality value'
            call cpu_time(end_time)
            result%elapsed = end_time - start_time
            return
         end if
      end if
      y = max(lower, min(upper, y))
      lambda = 0.0_dp
      h_inv = eye(n + m)
      rho = max(ctrl%rho, 1.0e-6_dp)
      if (nc == 0) rho = 0.0_dp
      result%n_eval = 0
      total_inner = 0

      if (nc > 0 .and. ctrl%restoration_iter > 0) then
         call restore_feasibility(prob, y, lower, upper, ctrl, result%n_eval)
      end if

      call raw_state(prob, y, f, c, ok, result%n_eval)
      if (.not. ok) then
         result%convergence = solnp_numerical_failure
         result%message = 'nonfinite objective or constraints at initial point'
         call cpu_time(end_time)
         result%elapsed = end_time - start_time
         return
      end if
      previous_f = f
      previous_cnorm = vnorm2(c)
      obj_hist(1) = f
      con_hist(1) = previous_cnorm
      step_hist(1) = 0.0_dp

      do outer = 1, ctrl%max_iter
         call minimize_augmented(prob, y, lower, upper, lambda, rho, h_inv, ctrl, &
            inner_used, step_norm, ok, result%n_eval)
         total_inner = total_inner + inner_used
         if (.not. ok) then
            result%convergence = solnp_numerical_failure
            result%message = 'inner augmented-Lagrangian minimization failed'
            exit
         end if

         call raw_state(prob, y, f, c, ok, result%n_eval)
         if (.not. ok) then
            result%convergence = solnp_numerical_failure
            result%message = 'nonfinite objective or constraints during outer iteration'
            exit
         end if
         cnorm = vnorm2(c)
         rel_obj = abs(f - previous_f) / max(1.0_dp, abs(previous_f))
         obj_hist(outer + 1) = f
         con_hist(outer + 1) = cnorm
         step_hist(outer + 1) = step_norm

         if (ctrl%trace > 0) then
            write(*, '(a,i0,2(a,es13.5),a,es10.3,a,es10.3)') &
               'iter ', outer, ' objective ', f, ' constraint ', cnorm, &
               ' rel_obj ', rel_obj, ' rho ', rho
         end if

         if (nc > 0) then
            lambda = lambda - rho * c
            call estimate_lagrange(prob, y, ctrl%delta, ctrl%tol, lambda, ok, result%n_eval)
            if (.not. ok) lambda = lambda - rho * c
         end if
         call kkt_from_state(prob, y, lambda, ctrl%delta, ctrl%tol, result%kkt)
         if (cnorm <= max(ctrl%tol, 10.0_dp * epsilon(1.0_dp)) .and. &
             result%kkt%stationarity <= max(sqrt(ctrl%tol), 1.0e-6_dp) .and. &
             (rel_obj <= sqrt(ctrl%tol) .or. step_norm <= sqrt(ctrl%tol))) then
            result%convergence = solnp_success
            result%message = 'converged'
            exit
         end if
         if (nc == 0 .and. result%kkt%stationarity <= max(sqrt(ctrl%tol), 1.0e-6_dp)) then
            result%convergence = solnp_success
            result%message = 'converged'
            exit
         end if

         if (nc > 0) then
            if (cnorm > 0.5_dp * previous_cnorm) then
               rho = min(ctrl%max_rho, rho * ctrl%penalty_growth)
               h_inv = eye(n + m)
            end if
         end if
         previous_cnorm = max(cnorm, tiny(1.0_dp))
         previous_f = f
         if (outer == ctrl%max_iter) then
            if (cnorm > max(100.0_dp * ctrl%tol, 1.0e-5_dp)) then
               result%convergence = solnp_infeasible
               result%message = 'maximum iterations reached with constraint violation'
            else
               result%convergence = solnp_max_iterations
               result%message = 'maximum outer iterations reached'
            end if
         end if
      end do

      if (result%convergence == solnp_invalid_problem) then
         result%convergence = solnp_max_iterations
         result%message = 'maximum outer iterations reached'
      end if
      result%out_iterations = min(outer, ctrl%max_iter)
      result%inner_iterations = total_inner
      allocate(result%pars(n), result%lagrange(nc), result%ineq_slack(m))
      result%pars = y(m + 1:)
      if (nc > 0) result%lagrange = lambda
      if (m > 0) result%ineq_slack = y(1:m)
      call evaluate_objective(prob, result%pars, result%objective, ok)
      result%n_eval = result%n_eval + 1
      allocate(result%objective_history(result%out_iterations + 1))
      allocate(result%constraint_history(result%out_iterations + 1))
      allocate(result%step_history(result%out_iterations + 1))
      result%objective_history = obj_hist(1:result%out_iterations + 1)
      result%constraint_history = con_hist(1:result%out_iterations + 1)
      result%step_history = step_hist(1:result%out_iterations + 1)

      allocate(full_hess(n + m, n + m), result%hessian(n, n))
      call invert_matrix(h_inv, full_hess, inv_ok)
      if (inv_ok) then
         result%hessian = full_hess(m + 1:, m + 1:)
      else
         result%hessian = eye(n)
      end if
      call symmetrize(result%hessian)
      call kkt_diagnose(prob, result, result%kkt, ctrl%tol)
      call cpu_time(end_time)
      result%elapsed = end_time - start_time
   end subroutine csolnp

   subroutine initialize_slack(problem, x, slack, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: slack(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: h(:)

      allocate(h(problem%n_ineq))
      call evaluate_inequalities(problem, x, h, ok)
      if (ok) slack = max(problem%ineq_lower, min(problem%ineq_upper, h))
   end subroutine initialize_slack

   subroutine raw_state(problem, y, f, c, ok, n_eval)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: c(:)
      logical, intent(out) :: ok
      integer, intent(inout) :: n_eval
      integer :: m

      m = problem%n_ineq
      call evaluate_objective(problem, y(m + 1:), f, ok)
      n_eval = n_eval + 1
      if (.not. ok) return
      call constraint_values(problem, y(m + 1:), y(1:m), c, ok)
   end subroutine raw_state

   subroutine augmented_value_gradient(problem, y, lambda, rho, value, gradient, c, ok, &
                                       delta, n_eval)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: y(:), lambda(:), rho, delta
      real(dp), intent(out) :: value, gradient(:), c(:)
      logical, intent(out) :: ok
      integer, intent(inout) :: n_eval

      real(dp), allocatable :: gf(:), j(:, :), q(:)
      real(dp) :: f
      logical :: ok1, ok2
      integer :: n, m, nc

      n = problem%n
      m = problem%n_ineq
      nc = problem%n_eq + problem%n_ineq
      allocate(gf(n), j(nc, n), q(nc))
      call evaluate_objective(problem, y(m + 1:), f, ok1)
      n_eval = n_eval + 1
      if (.not. associated(problem%gr)) n_eval = n_eval + 2 * n
      call evaluate_gradient(problem, y(m + 1:), gf, delta, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) return
      call constraint_values(problem, y(m + 1:), y(1:m), c, ok1)
      call constraint_jacobian(problem, y(m + 1:), j, delta, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) return

      value = f
      gradient = 0.0_dp
      gradient(m + 1:) = gf
      if (nc > 0) then
         q = -lambda + rho * c
         value = value - dot_product(lambda, c) + 0.5_dp * rho * dot_product(c, c)
         if (m > 0) gradient(1:m) = -q(problem%n_eq + 1:)
         gradient(m + 1:) = gradient(m + 1:) + matmul(transpose(j), q)
      end if
      ok = ieee_is_finite(value) .and. all(ieee_is_finite(gradient))
   end subroutine augmented_value_gradient

   subroutine augmented_value(problem, y, lambda, rho, value, ok, n_eval)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: y(:), lambda(:), rho
      real(dp), intent(out) :: value
      logical, intent(out) :: ok
      integer, intent(inout) :: n_eval
      real(dp), allocatable :: c(:)
      real(dp) :: f
      integer :: m, nc

      m = problem%n_ineq
      nc = problem%n_eq + problem%n_ineq
      allocate(c(nc))
      call evaluate_objective(problem, y(m + 1:), f, ok)
      n_eval = n_eval + 1
      if (.not. ok) return
      call constraint_values(problem, y(m + 1:), y(1:m), c, ok)
      if (.not. ok) return
      value = f
      if (nc > 0) value = value - dot_product(lambda, c) + 0.5_dp * rho * dot_product(c, c)
      ok = ieee_is_finite(value)
   end subroutine augmented_value

   subroutine minimize_augmented(problem, y, lower, upper, lambda, rho, h_inv, ctrl, &
                                 iterations, step_norm, ok, n_eval)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(inout) :: y(:), h_inv(:, :)
      real(dp), intent(in) :: lower(:), upper(:), lambda(:), rho
      type(solnp_control), intent(in) :: ctrl
      integer, intent(out) :: iterations
      real(dp), intent(out) :: step_norm
      logical, intent(out) :: ok
      integer, intent(inout) :: n_eval

      real(dp), allocatable :: g(:), gnew(:), pg(:), direction(:), trial(:), step(:)
      real(dp), allocatable :: c(:), cnew(:), ydiff(:), v(:, :)
      real(dp) :: value, trial_value, alpha, slope, ys, inv_ys, direction_norm, max_direction
      logical :: eval_ok, accepted
      integer :: k, ls, dim, nc

      dim = size(y)
      nc = problem%n_eq + problem%n_ineq
      allocate(g(dim), gnew(dim), pg(dim), direction(dim), trial(dim), step(dim))
      allocate(c(nc), cnew(nc), ydiff(dim), v(dim, dim))
      call augmented_value_gradient(problem, y, lambda, rho, value, g, c, ok, ctrl%delta, n_eval)
      if (.not. ok) then
         iterations = 0
         step_norm = 0.0_dp
         return
      end if
      step_norm = 0.0_dp
      do k = 1, ctrl%min_iter
         call projected_gradient(y, g, lower, upper, pg)
         if (norm_inf(pg) <= max(0.1_dp * sqrt(ctrl%tol), 1.0e-7_dp)) exit
         direction = -matmul(h_inv, pg)
         if (dot_product(direction, pg) >= -1.0e-12_dp * max(1.0_dp, vnorm2(direction) * vnorm2(pg))) then
            direction = -pg
            h_inv = eye(dim)
         end if
         direction_norm = vnorm2(direction)
         max_direction = 10.0_dp * max(1.0_dp, vnorm2(y))
         if (direction_norm > max_direction) direction = direction * (max_direction / direction_norm)
         alpha = 1.0_dp
         accepted = .false.
         do ls = 1, ctrl%line_search_max
            trial = max(lower, min(upper, y + alpha * direction))
            step = trial - y
            if (vnorm2(step) <= ctrl%min_step) then
               alpha = 0.5_dp * alpha
               cycle
            end if
            slope = dot_product(g, step)
            if (slope >= 0.0_dp) then
               direction = -pg
               trial = max(lower, min(upper, y + alpha * direction))
               step = trial - y
               slope = dot_product(g, step)
            end if
            call augmented_value(problem, trial, lambda, rho, trial_value, eval_ok, n_eval)
            if (eval_ok) then
               if (trial_value <= value + ctrl%armijo * slope) then
                  accepted = .true.
                  exit
               end if
            end if
            alpha = 0.5_dp * alpha
         end do
         if (.not. accepted) then
            step = -max(ctrl%min_step, min(1.0e-4_dp, 1.0_dp / max(1.0_dp, vnorm2(pg)))) * pg
            trial = max(lower, min(upper, y + step))
            step = trial - y
            call augmented_value(problem, trial, lambda, rho, trial_value, eval_ok, n_eval)
            if (.not. eval_ok .or. trial_value >= value) then
               exit
            end if
         end if
         call augmented_value_gradient(problem, trial, lambda, rho, trial_value, gnew, cnew, &
            eval_ok, ctrl%delta, n_eval)
         if (.not. eval_ok) then
            ok = .false.
            iterations = k
            return
         end if
         step = trial - y
         ydiff = gnew - g
         ys = dot_product(ydiff, step)
         if (ys > 1.0e-12_dp * max(1.0_dp, vnorm2(step) * vnorm2(ydiff))) then
            inv_ys = 1.0_dp / ys
            v = eye(dim) - inv_ys * outer_product(step, ydiff)
            h_inv = matmul(v, matmul(h_inv, transpose(v))) + inv_ys * outer_product(step, step)
            call symmetrize(h_inv)
         else if (mod(k, 20) == 0) then
            h_inv = eye(dim)
         end if
         y = trial
         g = gnew
         c = cnew
         value = trial_value
         step_norm = vnorm2(step)
         if (step_norm <= max(ctrl%min_step, 0.1_dp * ctrl%tol * max(1.0_dp, vnorm2(y)))) exit
      end do
      iterations = min(k, ctrl%min_iter)
      ok = .true.
   end subroutine minimize_augmented

   subroutine restore_feasibility(problem, y, lower, upper, ctrl, n_eval)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(inout) :: y(:)
      real(dp), intent(in) :: lower(:), upper(:)
      type(solnp_control), intent(in) :: ctrl
      integer, intent(inout) :: n_eval

      real(dp), allocatable :: c(:), ctrial(:), j(:, :), g(:), pg(:), trial(:)
      real(dp) :: value, trial_value, alpha
      logical :: ok1, ok2
      integer :: k, ls, n, m, nc

      n = problem%n
      m = problem%n_ineq
      nc = problem%n_eq + problem%n_ineq
      allocate(c(nc), ctrial(nc), j(nc, n), g(n + m), pg(n + m), trial(n + m))
      do k = 1, ctrl%restoration_iter
         call constraint_values(problem, y(m + 1:), y(1:m), c, ok1)
         call constraint_jacobian(problem, y(m + 1:), j, ctrl%delta, ok2)
         if (.not. ok1 .or. .not. ok2) return
         value = 0.5_dp * dot_product(c, c)
         if (sqrt(2.0_dp * value) <= max(ctrl%tol, 1.0e-8_dp)) return
         g = 0.0_dp
         if (m > 0) g(1:m) = -c(problem%n_eq + 1:)
         g(m + 1:) = matmul(transpose(j), c)
         call projected_gradient(y, g, lower, upper, pg)
         if (norm_inf(pg) <= 1.0e-10_dp) return
         alpha = min(1.0_dp, 1.0_dp / max(1.0_dp, vnorm2(pg)))
         do ls = 1, 24
            trial = max(lower, min(upper, y - alpha * pg))
            call constraint_values(problem, trial(m + 1:), trial(1:m), ctrial, ok1)
            if (ok1) then
               trial_value = 0.5_dp * dot_product(ctrial, ctrial)
               if (trial_value <= value - 1.0e-4_dp * alpha * dot_product(pg, pg)) exit
            end if
            alpha = 0.5_dp * alpha
         end do
         if (ls > 24) return
         y = trial
      end do
      n_eval = n_eval
   end subroutine restore_feasibility

   subroutine estimate_lagrange(problem, y, delta, tol, lambda, ok, n_eval)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: y(:), delta, tol
      real(dp), intent(inout) :: lambda(:)
      logical, intent(out) :: ok
      integer, intent(inout) :: n_eval

      real(dp), allocatable :: gradient(:), jacobian(:, :), active_j(:, :), normal(:, :)
      real(dp), allocatable :: normal_inv(:, :), rhs(:), solution(:)
      integer, allocatable :: active_index(:), active_side(:)
      real(dp) :: active_tol, regularization
      integer :: n, m, nc, na, i, row
      logical :: grad_ok, jac_ok, inverse_ok

      n = problem%n
      m = problem%n_ineq
      nc = problem%n_eq + m
      ok = .false.
      if (nc == 0) then
         ok = .true.
         return
      end if
      allocate(gradient(n), jacobian(nc, n), active_index(nc), active_side(nc))
      call evaluate_gradient(problem, y(m + 1:), gradient, delta, grad_ok)
      if (.not. associated(problem%gr)) n_eval = n_eval + 2 * n
      call constraint_jacobian(problem, y(m + 1:), jacobian, delta, jac_ok)
      if (.not. grad_ok .or. .not. jac_ok) return

      na = problem%n_eq
      do i = 1, problem%n_eq
         active_index(i) = i
         active_side(i) = 0
      end do
      active_tol = max(sqrt(tol), 1.0e-6_dp)
      do i = 1, m
         row = problem%n_eq + i
         if (y(i) <= problem%ineq_lower(i) + active_tol * &
             max(1.0_dp, abs(problem%ineq_lower(i)))) then
            na = na + 1
            active_index(na) = row
            active_side(na) = 1
         else if (y(i) >= problem%ineq_upper(i) - active_tol * &
                  max(1.0_dp, abs(problem%ineq_upper(i)))) then
            na = na + 1
            active_index(na) = row
            active_side(na) = -1
         end if
      end do
      lambda = 0.0_dp
      if (na == 0) then
         ok = .true.
         return
      end if
      allocate(active_j(na, n), normal(na, na), normal_inv(na, na), rhs(na), solution(na))
      do i = 1, na
         active_j(i, :) = jacobian(active_index(i), :)
      end do
      normal = matmul(active_j, transpose(active_j))
      regularization = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(normal)))
      do i = 1, na
         normal(i, i) = normal(i, i) + regularization
      end do
      rhs = matmul(active_j, gradient)
      call invert_matrix(normal, normal_inv, inverse_ok)
      if (.not. inverse_ok) return
      solution = matmul(normal_inv, rhs)
      do i = 1, na
         row = active_index(i)
         select case (active_side(i))
         case (1)
            lambda(row) = max(0.0_dp, solution(i))
         case (-1)
            lambda(row) = min(0.0_dp, solution(i))
         case default
            lambda(row) = solution(i)
         end select
      end do
      ok = all(ieee_is_finite(lambda))
   end subroutine estimate_lagrange

   subroutine kkt_diagnose(problem, result, diagnostics, tol)
      type(solnp_problem), intent(in) :: problem
      type(solnp_result), intent(in) :: result
      type(kkt_diagnostics), intent(out) :: diagnostics
      real(dp), intent(in), optional :: tol

      type(solnp_problem) :: prob
      real(dp), allocatable :: y(:), lambda(:)
      real(dp) :: tolerance
      integer :: stat, m, nc
      character(len=160) :: msg

      tolerance = 1.0e-6_dp
      if (present(tol)) tolerance = tol
      prob = problem
      call prepare_problem(prob, stat, msg)
      if (stat /= solnp_success .or. .not. allocated(result%pars)) then
         diagnostics = kkt_diagnostics()
         return
      end if
      m = prob%n_ineq
      nc = prob%n_eq + prob%n_ineq
      allocate(y(prob%n + m), lambda(nc))
      y(m + 1:) = result%pars
      if (m > 0) then
         if (allocated(result%ineq_slack)) then
            y(1:m) = result%ineq_slack
         else
            call initialize_slack(prob, result%pars, y(1:m), diagnostics%primal_feasible)
         end if
      end if
      lambda = 0.0_dp
      if (allocated(result%lagrange)) then
         if (size(result%lagrange) == nc) lambda = result%lagrange
      end if
      call kkt_from_state(prob, y, lambda, max(1.0e-7_dp, sqrt(epsilon(1.0_dp))), &
         tolerance, diagnostics)
   end subroutine kkt_diagnose

   subroutine kkt_from_state(problem, y, lambda, delta, tol, diagnostics)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: y(:), lambda(:), delta, tol
      type(kkt_diagnostics), intent(out) :: diagnostics

      real(dp), allocatable :: gf(:), eq(:), h(:), jeq(:, :), jineq(:, :), lag_grad(:), pg(:)
      real(dp) :: lower_gap, upper_gap, dual_v, comp
      logical :: ok1, ok2
      integer :: n, m, i

      n = problem%n
      m = problem%n_ineq
      allocate(gf(n), eq(problem%n_eq), h(m), jeq(problem%n_eq, n), jineq(m, n))
      allocate(lag_grad(n), pg(n))
      call evaluate_gradient(problem, y(m + 1:), gf, delta, ok1)
      call evaluate_equalities(problem, y(m + 1:), eq, ok2)
      if (.not. ok1 .or. .not. ok2) then
         diagnostics = kkt_diagnostics()
         return
      end if
      call evaluate_inequalities(problem, y(m + 1:), h, ok1)
      call evaluate_eq_jacobian(problem, y(m + 1:), jeq, delta, ok2)
      if (.not. ok1 .or. .not. ok2) then
         diagnostics = kkt_diagnostics()
         return
      end if
      call evaluate_ineq_jacobian(problem, y(m + 1:), jineq, delta, ok1)
      if (.not. ok1) then
         diagnostics = kkt_diagnostics()
         return
      end if
      lag_grad = gf
      if (problem%n_eq > 0) lag_grad = lag_grad - &
         matmul(transpose(jeq), lambda(1:problem%n_eq))
      if (m > 0) lag_grad = lag_grad - &
         matmul(transpose(jineq), lambda(problem%n_eq + 1:))
      call projected_gradient(y(m + 1:), lag_grad, problem%lower, problem%upper, pg)
      diagnostics%stationarity = norm_inf(pg)
      diagnostics%eq_violation = norm_inf(eq)
      diagnostics%ineq_violation = 0.0_dp
      diagnostics%bound_violation = max(maxval(max(problem%lower - y(m + 1:), 0.0_dp)), &
         maxval(max(y(m + 1:) - problem%upper, 0.0_dp)))
      dual_v = 0.0_dp
      comp = 0.0_dp
      do i = 1, m
         diagnostics%ineq_violation = max(diagnostics%ineq_violation, &
            max(problem%ineq_lower(i) - h(i), h(i) - problem%ineq_upper(i), 0.0_dp))
         lower_gap = max(0.0_dp, y(i) - problem%ineq_lower(i))
         upper_gap = max(0.0_dp, problem%ineq_upper(i) - y(i))
         if (y(i) <= problem%ineq_lower(i) + sqrt(tol) * max(1.0_dp, abs(problem%ineq_lower(i)))) then
            dual_v = max(dual_v, max(0.0_dp, -lambda(problem%n_eq + i)))
            comp = max(comp, abs(lambda(problem%n_eq + i)) * lower_gap)
         else if (y(i) >= problem%ineq_upper(i) - sqrt(tol) * max(1.0_dp, abs(problem%ineq_upper(i)))) then
            dual_v = max(dual_v, max(0.0_dp, lambda(problem%n_eq + i)))
            comp = max(comp, abs(lambda(problem%n_eq + i)) * upper_gap)
         else
            dual_v = max(dual_v, abs(lambda(problem%n_eq + i)))
            comp = max(comp, abs(lambda(problem%n_eq + i)) * min(lower_gap, upper_gap))
         end if
      end do
      diagnostics%dual_feas_violation = dual_v
      diagnostics%complementarity = comp
      diagnostics%primal_feasible = max(diagnostics%eq_violation, diagnostics%ineq_violation, &
         diagnostics%bound_violation) <= max(tol, 1.0e-6_dp)
      diagnostics%first_order = diagnostics%primal_feasible .and. &
         diagnostics%stationarity <= max(sqrt(tol), 1.0e-5_dp) .and. &
         diagnostics%dual_feas_violation <= max(sqrt(tol), 1.0e-5_dp)
   end subroutine kkt_from_state

end module rsolnp_solver
