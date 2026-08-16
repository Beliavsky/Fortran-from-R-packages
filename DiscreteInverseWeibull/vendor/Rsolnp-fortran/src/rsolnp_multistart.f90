! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_multistart
   use rsolnp_kinds, only : dp
   use rsolnp_types, only : solnp_problem, solnp_control, solnp_result, multistart_result, &
      solnp_success
   use rsolnp_problem, only : prepare_problem
   use rsolnp_solver, only : csolnp
   use rsolnp_evaluate, only : evaluate_objective, evaluate_equalities, evaluate_inequalities, &
      evaluate_eq_jacobian, evaluate_ineq_jacobian
   use rsolnp_linalg, only : vnorm2, norm_inf, projected_gradient
   implicit none
   private

   public :: startpars, csolnp_ms, gosolnp

contains

   subroutine startpars(problem, n_starts, starts, seed, include_start, feasibility_iter, status, message)
      type(solnp_problem), intent(in) :: problem
      integer, intent(in) :: n_starts
      real(dp), allocatable, intent(out) :: starts(:, :)
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: include_start
      integer, intent(in), optional :: feasibility_iter
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message

      type(solnp_problem) :: prob
      real(dp), allocatable :: lo(:), hi(:), merit(:), temp(:)
      real(dp) :: u
      integer :: i, j, offset, stat, feas_iter, best
      logical :: keep_start
      character(len=160) :: msg

      prob = problem
      call prepare_problem(prob, stat, msg)
      if (stat /= solnp_success .or. n_starts < 1) then
         allocate(starts(0, 0))
         if (present(status)) status = merge(stat, 2, stat /= solnp_success)
         if (present(message)) then
            if (stat /= solnp_success) then
               message = trim(msg)
            else
               message = 'n_starts must be positive'
            end if
         end if
         return
      end if
      keep_start = .true.
      if (present(include_start)) keep_start = include_start
      feas_iter = 30
      if (present(feasibility_iter)) feas_iter = max(0, feasibility_iter)
      offset = 0
      if (present(seed)) offset = abs(seed)
      allocate(starts(n_starts, prob%n), lo(prob%n), hi(prob%n), merit(n_starts), temp(prob%n))
      call practical_bounds(prob, lo, hi)
      do i = 1, n_starts
         do j = 1, prob%n
            if (abs(prob%lower(j) - prob%upper(j)) <= epsilon(1.0_dp) * max(1.0_dp, abs(prob%lower(j)))) then
               starts(i, j) = prob%lower(j)
            else
               u = halton(i + offset, prime_number(j))
               starts(i, j) = lo(j) + u * (hi(j) - lo(j))
            end if
         end do
      end do
      if (keep_start) starts(1, :) = prob%start
      if (feas_iter > 0 .and. (prob%n_eq + prob%n_ineq > 0)) then
         do i = 1, n_starts
            temp = starts(i, :)
            call quick_feasible(prob, temp, feas_iter)
            starts(i, :) = temp
         end do
      end if
      do i = 1, n_starts
         merit(i) = start_merit(prob, starts(i, :))
      end do
      do i = 1, n_starts - 1
         best = i - 1 + minloc(merit(i:), dim=1)
         if (best /= i) then
            u = merit(i)
            merit(i) = merit(best)
            merit(best) = u
            temp = starts(i, :)
            starts(i, :) = starts(best, :)
            starts(best, :) = temp
         end if
      end do
      if (present(status)) status = solnp_success
      if (present(message)) message = 'success'
   end subroutine startpars

   subroutine csolnp_ms(problem, n_starts, result, control, seed)
      type(solnp_problem), intent(in) :: problem
      integer, intent(in) :: n_starts
      type(multistart_result), intent(out) :: result
      type(solnp_control), intent(in), optional :: control
      integer, intent(in), optional :: seed

      type(solnp_problem) :: prob_i, prob
      type(solnp_control) :: ctrl
      integer :: i, stat, best_index
      real(dp) :: score, best_score, violation
      character(len=160) :: msg

      ctrl = solnp_control()
      if (present(control)) ctrl = control
      prob = problem
      call prepare_problem(prob, stat, msg)
      if (stat /= solnp_success .or. n_starts < 1) then
         result%best%convergence = 2
         result%best%message = trim(msg)
         return
      end if
      call startpars(prob, n_starts, result%starts, seed=seed, status=stat, message=msg)
      if (stat /= solnp_success) then
         result%best%convergence = stat
         result%best%message = trim(msg)
         return
      end if
      allocate(result%results(n_starts), result%objectives(n_starts), result%convergence(n_starts))
      best_score = huge(1.0_dp)
      best_index = 0
      do i = 1, n_starts
         prob_i = prob
         prob_i%start = result%starts(i, :)
         call csolnp(prob_i, result%results(i), ctrl)
         result%objectives(i) = result%results(i)%objective
         result%convergence(i) = result%results(i)%convergence
         violation = max(result%results(i)%kkt%eq_violation, &
            result%results(i)%kkt%ineq_violation, result%results(i)%kkt%bound_violation)
         score = result%results(i)%objective + 1.0e8_dp * violation
         if (result%results(i)%convergence /= solnp_success) score = score + 1.0e6_dp
         if (score < best_score) then
            best_score = score
            best_index = i
         end if
      end do
      result%best_index = best_index
      if (best_index > 0) result%best = result%results(best_index)
   end subroutine csolnp_ms

   subroutine gosolnp(problem, n_starts, result, control, seed)
      type(solnp_problem), intent(in) :: problem
      integer, intent(in) :: n_starts
      type(multistart_result), intent(out) :: result
      type(solnp_control), intent(in), optional :: control
      integer, intent(in), optional :: seed
      call csolnp_ms(problem, n_starts, result, control, seed)
   end subroutine gosolnp

   subroutine practical_bounds(problem, lower, upper)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(out) :: lower(:), upper(:)
      real(dp) :: scale
      integer :: i

      do i = 1, problem%n
         scale = 10.0_dp * max(1.0_dp, abs(problem%start(i)))
         lower(i) = problem%lower(i)
         upper(i) = problem%upper(i)
         if (lower(i) < -1.0e12_dp) lower(i) = problem%start(i) - scale
         if (upper(i) > 1.0e12_dp) upper(i) = problem%start(i) + scale
         if (abs(lower(i) - upper(i)) <= epsilon(1.0_dp) * max(1.0_dp, abs(lower(i)))) cycle
         if (upper(i) - lower(i) > 1.0e6_dp * max(1.0_dp, scale)) then
            lower(i) = max(lower(i), problem%start(i) - scale)
            upper(i) = min(upper(i), problem%start(i) + scale)
         end if
      end do
   end subroutine practical_bounds

   subroutine quick_feasible(problem, x, max_iter)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(inout) :: x(:)
      integer, intent(in) :: max_iter

      real(dp), allocatable :: eq(:), h(:), jeq(:, :), jineq(:, :), residual(:), gradient(:), pg(:), trial(:)
      real(dp) :: value, trial_value, alpha
      logical :: ok1, ok2
      integer :: k, ls, nc

      nc = problem%n_eq + problem%n_ineq
      allocate(eq(problem%n_eq), h(problem%n_ineq), jeq(problem%n_eq, problem%n))
      allocate(jineq(problem%n_ineq, problem%n), residual(nc), gradient(problem%n), &
         pg(problem%n), trial(problem%n))
      do k = 1, max_iter
         call feasibility_state(problem, x, eq, h, residual, jeq, jineq, ok1)
         if (.not. ok1) return
         value = 0.5_dp * dot_product(residual, residual)
         if (sqrt(2.0_dp * value) <= 1.0e-8_dp) return
         gradient = 0.0_dp
         if (problem%n_eq > 0) gradient = gradient + matmul(transpose(jeq), residual(1:problem%n_eq))
         if (problem%n_ineq > 0) gradient = gradient + &
            matmul(transpose(jineq), residual(problem%n_eq + 1:))
         call projected_gradient(x, gradient, problem%lower, problem%upper, pg)
         if (norm_inf(pg) <= 1.0e-10_dp) return
         alpha = min(1.0_dp, 1.0_dp / max(1.0_dp, vnorm2(pg)))
         do ls = 1, 24
            trial = max(problem%lower, min(problem%upper, x - alpha * pg))
            call feasibility_residual(problem, trial, residual, ok2)
            if (ok2) then
               trial_value = 0.5_dp * dot_product(residual, residual)
               if (trial_value < value) exit
            end if
            alpha = 0.5_dp * alpha
         end do
         if (ls > 24) return
         x = trial
      end do
   end subroutine quick_feasible

   subroutine feasibility_state(problem, x, eq, h, residual, jeq, jineq, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: eq(:), h(:), residual(:), jeq(:, :), jineq(:, :)
      logical, intent(out) :: ok
      logical :: ok1, ok2
      integer :: i

      call evaluate_equalities(problem, x, eq, ok1)
      call evaluate_inequalities(problem, x, h, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) return
      call evaluate_eq_jacobian(problem, x, jeq, 1.0e-6_dp, ok1)
      call evaluate_ineq_jacobian(problem, x, jineq, 1.0e-6_dp, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) return
      if (problem%n_eq > 0) residual(1:problem%n_eq) = eq
      do i = 1, problem%n_ineq
         if (h(i) < problem%ineq_lower(i)) then
            residual(problem%n_eq + i) = h(i) - problem%ineq_lower(i)
         else if (h(i) > problem%ineq_upper(i)) then
            residual(problem%n_eq + i) = h(i) - problem%ineq_upper(i)
         else
            residual(problem%n_eq + i) = 0.0_dp
            jineq(i, :) = 0.0_dp
         end if
      end do
   end subroutine feasibility_state

   subroutine feasibility_residual(problem, x, residual, ok)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: residual(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: eq(:), h(:)
      logical :: ok1, ok2
      integer :: i

      allocate(eq(problem%n_eq), h(problem%n_ineq))
      call evaluate_equalities(problem, x, eq, ok1)
      call evaluate_inequalities(problem, x, h, ok2)
      ok = ok1 .and. ok2
      if (.not. ok) return
      if (problem%n_eq > 0) residual(1:problem%n_eq) = eq
      do i = 1, problem%n_ineq
         if (h(i) < problem%ineq_lower(i)) then
            residual(problem%n_eq + i) = h(i) - problem%ineq_lower(i)
         else if (h(i) > problem%ineq_upper(i)) then
            residual(problem%n_eq + i) = h(i) - problem%ineq_upper(i)
         else
            residual(problem%n_eq + i) = 0.0_dp
         end if
      end do
   end subroutine feasibility_residual

   real(dp) function start_merit(problem, x) result(value)
      type(solnp_problem), intent(in) :: problem
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: f
      logical :: ok1, ok2

      allocate(residual(problem%n_eq + problem%n_ineq))
      call evaluate_objective(problem, x, f, ok1)
      call feasibility_residual(problem, x, residual, ok2)
      if (.not. ok1 .or. .not. ok2) then
         value = huge(1.0_dp)
      else
         value = f + 1.0e6_dp * dot_product(residual, residual)
      end if
   end function start_merit

   integer function prime_number(index) result(p)
      integer, intent(in) :: index
      integer :: candidate, count, d
      logical :: is_prime

      count = 0
      candidate = 1
      do while (count < index)
         candidate = candidate + 1
         is_prime = .true.
         do d = 2, int(sqrt(real(candidate, dp)))
            if (mod(candidate, d) == 0) then
               is_prime = .false.
               exit
            end if
         end do
         if (is_prime) count = count + 1
      end do
      p = candidate
   end function prime_number

   pure real(dp) function halton(index, base) result(value)
      integer, intent(in) :: index, base
      integer :: i
      real(dp) :: factor

      i = index
      factor = 1.0_dp
      value = 0.0_dp
      do while (i > 0)
         factor = factor / real(base, dp)
         value = value + factor * real(mod(i, base), dp)
         i = i / base
      end do
   end function halton

end module rsolnp_multistart
