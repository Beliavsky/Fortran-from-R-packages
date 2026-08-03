! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_solvers
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nloptr_kinds, only: dp
  use nloptr_types
  use nloptr_utils, only: project_bounds, max_abs, vec_norm, outer_product, halton_point, &
    status_message, validate_problem
  use nloptr_evaluation, only: evaluate_penalty, evaluate_problem
  implicit none
  private
  public :: optimize_problem

contains

  pure function upper_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c
    out = text
    do i = 1, len(text)
      c = iachar(text(i:i))
      if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - 32)
    end do
  end function upper_string

  pure logical function is_global_algorithm(algorithm) result(answer)
    character(len=*), intent(in) :: algorithm
    character(len=len(algorithm)) :: alg
    alg = upper_string(algorithm)
    answer = index(alg, 'DIRECT') > 0 .or. index(alg, 'CRS') > 0 .or. &
      index(alg, 'ISRES') > 0 .or. index(alg, 'STOGO') > 0 .or. index(alg, 'MLSL') > 0
  end function is_global_algorithm

  pure logical function is_derivative_free(algorithm) result(answer)
    character(len=*), intent(in) :: algorithm
    character(len=len(algorithm)) :: alg
    alg = upper_string(algorithm)
    answer = index(alg, '_LN_') > 0 .or. index(alg, '_GN_') > 0 .or. &
      index(alg, 'NELDERMEAD') > 0 .or. index(alg, 'SBPLX') > 0 .or. &
      index(alg, 'COBYLA') > 0 .or. index(alg, 'BOBYQA') > 0 .or. &
      index(alg, 'NEWUOA') > 0 .or. index(alg, 'PRAXIS') > 0
  end function is_derivative_free

  subroutine optimize_problem(problem, x0, options, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in) :: options
    type(nloptr_result), intent(out) :: result
    real(dp), allocatable :: x(:)
    real(dp) :: rho, violation, f, dummy_grad(problem%n)
    real(dp), allocatable :: g(:), h(:), jg(:, :), jh(:, :)
    integer :: outer, status, evals_before
    logical :: derivative_free

    allocate(result%solution(problem%n), x(problem%n))
    status = validate_problem(problem, x0)
    if (status /= 0) then
      result%solution = x0
      result%status = status
      result%message = status_message(status)
      return
    end if

    x = x0
    call project_bounds(x, problem%lower, problem%upper)
    rho = max(1.0e-6_dp, options%penalty_initial)
    derivative_free = is_derivative_free(options%algorithm)

    if (is_global_algorithm(options%algorithm)) then
      call global_search(problem, x, options, rho, result)
      return
    end if

    result%evaluations = 0
    result%iterations = 0
    status = NLOPT_MAXEVAL_REACHED
    do outer = 1, max(1, options%max_outer)
      evals_before = result%evaluations
      if (derivative_free) then
        if (index(upper_string(options%algorithm), 'NELDER') > 0 .or. &
            index(upper_string(options%algorithm), 'SBPLX') > 0 .or. &
            index(upper_string(options%algorithm), 'NEWUOA') > 0 .or. &
            index(upper_string(options%algorithm), 'BOBYQA') > 0) then
          call nelder_mead_core(problem, x, options, rho, result%evaluations, &
            result%iterations, status)
        else
          call pattern_search_core(problem, x, options, rho, result%evaluations, &
            result%iterations, status)
        end if
      else
        call bfgs_core(problem, x, options, rho, result%evaluations, result%iterations, status)
      end if

      call final_values(problem, x, f, violation, result%evaluations, status)
      if (status < 0) exit
      if (violation <= options%constraint_tol) then
        if (status == NLOPT_MAXEVAL_REACHED .and. &
            result%evaluations - evals_before < options%maxeval) status = NLOPT_XTOL_REACHED
        exit
      end if
      if (result%evaluations >= options%maxeval) then
        status = NLOPT_MAXEVAL_REACHED
        exit
      end if
      rho = rho * max(2.0_dp, options%penalty_growth)
    end do

    allocate(g(problem%n_ineq), h(problem%n_eq))
    allocate(jg(problem%n_ineq, problem%n), jh(problem%n_eq, problem%n))
    call evaluate_problem(problem, x, f, dummy_grad, g, jg, h, jh, .false., &
      outer, result%evaluations)
    if (outer /= 0) status = outer
    result%solution = x
    result%objective = f
    if (problem%n_ineq + problem%n_eq == 0) then
      result%max_constraint = 0.0_dp
    else
      result%max_constraint = 0.0_dp
      if (problem%n_ineq > 0) result%max_constraint = max(result%max_constraint, &
        maxval(max(0.0_dp, g)))
      if (problem%n_eq > 0) result%max_constraint = max(result%max_constraint, maxval(abs(h)))
    end if
    result%status = status
    result%converged = status > 0 .and. result%max_constraint <= &
      max(10.0_dp * options%constraint_tol, sqrt(epsilon(1.0_dp)))
    result%message = status_message(status)
  end subroutine optimize_problem

  subroutine final_values(problem, x, f, violation, evals, status)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, violation
    integer, intent(inout) :: evals
    integer, intent(out) :: status
    real(dp) :: grad(problem%n)
    real(dp), allocatable :: g(:), h(:), jg(:, :), jh(:, :)
    allocate(g(problem%n_ineq), h(problem%n_eq))
    allocate(jg(problem%n_ineq, problem%n), jh(problem%n_eq, problem%n))
    call evaluate_problem(problem, x, f, grad, g, jg, h, jh, .false., status, evals)
    violation = 0.0_dp
    if (status /= 0) return
    if (problem%n_ineq > 0) violation = max(violation, maxval(max(0.0_dp, g)))
    if (problem%n_eq > 0) violation = max(violation, maxval(abs(h)))
  end subroutine final_values

  subroutine bfgs_core(problem, x, options, rho, evals, iterations, status)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(inout) :: x(:)
    type(nloptr_options), intent(in) :: options
    real(dp), intent(in) :: rho
    integer, intent(inout) :: evals, iterations
    integer, intent(out) :: status
    integer :: n, i, stat, line_iter
    real(dp) :: f, f_new, f_old, violation, violation_new, alpha, dg, ys, relx, relf
    real(dp), allocatable :: grad(:), grad_new(:), hess_inv(:, :), direction(:)
    real(dp), allocatable :: x_new(:), s(:), y(:), identity(:, :), v(:, :)

    n = size(x)
    allocate(grad(n), grad_new(n), hess_inv(n, n), direction(n), x_new(n), s(n), y(n))
    allocate(identity(n, n), v(n, n))
    identity = 0.0_dp
    do i = 1, n
      identity(i, i) = 1.0_dp
    end do
    hess_inv = identity
    call evaluate_penalty(problem, x, rho, f, grad, violation, .true., stat, evals)
    if (stat /= 0) then
      status = stat
      return
    end if
    status = NLOPT_MAXEVAL_REACHED

    do while (evals < options%maxeval)
      iterations = iterations + 1
      if (vec_norm(grad) <= max(options%xtol_abs, 1.0e-10_dp)) then
        status = NLOPT_XTOL_REACHED
        return
      end if
      direction = -matmul(hess_inv, grad)
      if (dot_product(direction, grad) >= -1.0e-14_dp * max(1.0_dp, vec_norm(grad))) then
        direction = -grad
        hess_inv = identity
      end if
      alpha = 1.0_dp
      dg = dot_product(grad, direction)
      do line_iter = 1, 40
        x_new = x + alpha * direction
        call project_bounds(x_new, problem%lower, problem%upper)
        call evaluate_penalty(problem, x_new, rho, f_new, grad_new, violation_new, .true., stat, evals)
        if (stat == 0) then
          if (f_new <= f + 1.0e-4_dp * alpha * dg) exit
        end if
        alpha = 0.5_dp * alpha
      end do
      if (line_iter > 40 .or. alpha <= epsilon(1.0_dp)) then
        status = NLOPT_ROUNDOFF_LIMITED
        return
      end if
      f_old = f
      s = x_new - x
      y = grad_new - grad
      relx = max_abs(s) / max(1.0_dp, max_abs(x))
      relf = abs(f_new - f_old) / max(1.0_dp, abs(f_old))
      ys = dot_product(y, s)
      if (ys > 1.0e-12_dp * vec_norm(y) * max(vec_norm(s), 1.0e-30_dp)) then
        v = identity - outer_product(s, y) / ys
        hess_inv = matmul(v, matmul(hess_inv, transpose(v))) + outer_product(s, s) / ys
      else
        hess_inv = identity
      end if
      x = x_new
      grad = grad_new
      f = f_new
      violation = violation_new
      if (f <= options%stopval) then
        status = NLOPT_STOPVAL_REACHED
        return
      end if
      if (relf <= options%ftol_rel .or. abs(f_new - f_old) <= options%ftol_abs) then
        status = NLOPT_FTOL_REACHED
        return
      end if
      if (relx <= options%xtol_rel .or. max_abs(s) <= options%xtol_abs) then
        status = NLOPT_XTOL_REACHED
        return
      end if
    end do
  end subroutine bfgs_core

  subroutine nelder_mead_core(problem, x, options, rho, evals, iterations, status)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(inout) :: x(:)
    type(nloptr_options), intent(in) :: options
    real(dp), intent(in) :: rho
    integer, intent(inout) :: evals, iterations
    integer, intent(out) :: status
    integer :: n, j, stat
    real(dp) :: step, fr, fe, fc, vr, ve, vc, spread_x, spread_f
    real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:), xe(:), xc(:), grad(:)

    n = size(x)
    allocate(simplex(n, n + 1), values(n + 1), centroid(n), xr(n), xe(n), xc(n), grad(n))
    simplex(:, 1) = x
    do j = 1, n
      simplex(:, j + 1) = x
      step = options%initial_step
      if (step <= 0.0_dp) step = 0.05_dp * max(1.0_dp, abs(x(j)))
      if (ieee_is_finite(problem%lower(j))) then
        if (ieee_is_finite(problem%upper(j))) then
          if (abs(problem%lower(j)) < sqrt(huge(1.0_dp)) .and. &
              abs(problem%upper(j)) < sqrt(huge(1.0_dp))) then
            step = min(step, 0.1_dp * max(problem%upper(j) - problem%lower(j), options%xtol_abs))
          end if
        end if
      end if
      simplex(j, j + 1) = simplex(j, j + 1) + step
      call project_bounds(simplex(:, j + 1), problem%lower, problem%upper)
      if (max_abs(simplex(:, j + 1) - x) <= epsilon(1.0_dp)) then
        simplex(j, j + 1) = simplex(j, j + 1) - 2.0_dp * step
        call project_bounds(simplex(:, j + 1), problem%lower, problem%upper)
      end if
    end do
    do j = 1, n + 1
      call evaluate_penalty(problem, simplex(:, j), rho, values(j), grad, vr, .false., stat, evals)
      if (stat /= 0) then
        status = stat
        return
      end if
    end do

    status = NLOPT_MAXEVAL_REACHED
    do while (evals < options%maxeval)
      iterations = iterations + 1
      call sort_simplex(simplex, values)
      spread_f = maxval(abs(values - values(1))) / max(1.0_dp, abs(values(1)))
      spread_x = 0.0_dp
      do j = 2, n + 1
        spread_x = max(spread_x, max_abs(simplex(:, j) - simplex(:, 1)))
      end do
      if (spread_f <= options%ftol_rel) then
        status = NLOPT_FTOL_REACHED
        exit
      end if
      if (spread_x <= max(options%xtol_abs, options%xtol_rel * max(1.0_dp, max_abs(simplex(:, 1))))) then
        status = NLOPT_XTOL_REACHED
        exit
      end if
      centroid = sum(simplex(:, 1:n), dim=2) / real(n, dp)
      xr = centroid + (centroid - simplex(:, n + 1))
      call project_bounds(xr, problem%lower, problem%upper)
      call evaluate_penalty(problem, xr, rho, fr, grad, vr, .false., stat, evals)
      if (fr < values(1)) then
        xe = centroid + 2.0_dp * (xr - centroid)
        call project_bounds(xe, problem%lower, problem%upper)
        call evaluate_penalty(problem, xe, rho, fe, grad, ve, .false., stat, evals)
        if (fe < fr) then
          simplex(:, n + 1) = xe
          values(n + 1) = fe
        else
          simplex(:, n + 1) = xr
          values(n + 1) = fr
        end if
      else if (fr < values(n)) then
        simplex(:, n + 1) = xr
        values(n + 1) = fr
      else
        if (fr < values(n + 1)) then
          xc = centroid + 0.5_dp * (xr - centroid)
        else
          xc = centroid + 0.5_dp * (simplex(:, n + 1) - centroid)
        end if
        call project_bounds(xc, problem%lower, problem%upper)
        call evaluate_penalty(problem, xc, rho, fc, grad, vc, .false., stat, evals)
        if (fc < min(fr, values(n + 1))) then
          simplex(:, n + 1) = xc
          values(n + 1) = fc
        else
          do j = 2, n + 1
            simplex(:, j) = simplex(:, 1) + 0.5_dp * (simplex(:, j) - simplex(:, 1))
            call project_bounds(simplex(:, j), problem%lower, problem%upper)
            call evaluate_penalty(problem, simplex(:, j), rho, values(j), grad, vr, .false., stat, evals)
          end do
        end if
      end if
      if (values(1) <= options%stopval) then
        status = NLOPT_STOPVAL_REACHED
        exit
      end if
    end do
    call sort_simplex(simplex, values)
    x = simplex(:, 1)
  end subroutine nelder_mead_core

  subroutine sort_simplex(simplex, values)
    real(dp), intent(inout) :: simplex(:, :), values(:)
    integer :: i, j, best
    real(dp) :: temp
    real(dp) :: column(size(simplex, 1))
    do i = 1, size(values) - 1
      best = i
      do j = i + 1, size(values)
        if (values(j) < values(best)) best = j
      end do
      if (best /= i) then
        temp = values(i)
        values(i) = values(best)
        values(best) = temp
        column = simplex(:, i)
        simplex(:, i) = simplex(:, best)
        simplex(:, best) = column
      end if
    end do
  end subroutine sort_simplex

  subroutine pattern_search_core(problem, x, options, rho, evals, iterations, status)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(inout) :: x(:)
    type(nloptr_options), intent(in) :: options
    real(dp), intent(in) :: rho
    integer, intent(inout) :: evals, iterations
    integer, intent(out) :: status
    integer :: i, sign_direction, stat
    real(dp) :: f, trial_f, violation, trial_violation, step, scale
    real(dp) :: trial(size(x)), grad(size(x))
    logical :: improved

    scale = max(1.0_dp, max_abs(x))
    step = options%initial_step
    if (step <= 0.0_dp) step = 0.25_dp * scale
    call evaluate_penalty(problem, x, rho, f, grad, violation, .false., stat, evals)
    if (stat /= 0) then
      status = stat
      return
    end if
    status = NLOPT_MAXEVAL_REACHED
    do while (evals < options%maxeval)
      iterations = iterations + 1
      improved = .false.
      do i = 1, size(x)
        do sign_direction = -1, 1, 2
          trial = x
          trial(i) = trial(i) + real(sign_direction, dp) * step
          call project_bounds(trial, problem%lower, problem%upper)
          call evaluate_penalty(problem, trial, rho, trial_f, grad, trial_violation, .false., stat, evals)
          if (stat == 0 .and. trial_f < f) then
            x = trial
            f = trial_f
            violation = trial_violation
            improved = .true.
          end if
          if (evals >= options%maxeval) exit
        end do
        if (evals >= options%maxeval) exit
      end do
      if (.not. improved) step = 0.5_dp * step
      if (f <= options%stopval) then
        status = NLOPT_STOPVAL_REACHED
        return
      end if
      if (step <= max(options%xtol_abs, options%xtol_rel * max(1.0_dp, max_abs(x)))) then
        status = NLOPT_XTOL_REACHED
        return
      end if
    end do
  end subroutine pattern_search_core

  subroutine global_search(problem, x0, options, rho, result)
    type(nloptr_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(nloptr_options), intent(in) :: options
    real(dp), intent(in) :: rho
    type(nloptr_result), intent(out) :: result
    integer :: n, population, i, j, stat, local_evals, local_iters, n_local
    real(dp) :: value, violation, best_value, final_f, final_violation
    real(dp), allocatable :: points(:, :), values(:), unit(:), grad(:), candidate(:), best(:)
    type(nloptr_options) :: local_options

    n = size(x0)
    if (any(.not. ieee_is_finite(problem%lower)) .or. any(.not. ieee_is_finite(problem%upper))) then
      allocate(result%solution(n))
      result%solution = x0
      result%status = NLOPT_INVALID_ARGS
      result%message = 'global algorithms require finite lower and upper bounds'
      return
    end if
    population = options%population
    if (population <= 0) population = max(40, 12 * n)
    allocate(points(n, population + 1), values(population + 1), unit(n), grad(n), candidate(n), best(n))
    points(:, 1) = x0
    call project_bounds(points(:, 1), problem%lower, problem%upper)
    do i = 1, population
      call halton_point(i + abs(options%seed), unit)
      points(:, i + 1) = problem%lower + unit * (problem%upper - problem%lower)
    end do
    result%evaluations = 0
    result%iterations = 0
    do i = 1, population + 1
      call evaluate_penalty(problem, points(:, i), rho, values(i), grad, violation, .false., stat, result%evaluations)
      if (stat /= 0) values(i) = huge(1.0_dp)
    end do
    best_value = minval(values)
    best = points(:, minloc(values, dim=1))

    local_options = options
    local_options%algorithm = options%local_algorithm
    local_options%maxeval = max(100, (options%maxeval - result%evaluations) / 5)
    local_options%max_outer = 1
    n_local = min(5, population + 1)
    do j = 1, n_local
      i = minloc(values, dim=1)
      candidate = points(:, i)
      values(i) = huge(1.0_dp)
      local_evals = 0
      local_iters = 0
      if (is_derivative_free(local_options%algorithm)) then
        call nelder_mead_core(problem, candidate, local_options, rho, local_evals, local_iters, stat)
      else
        call bfgs_core(problem, candidate, local_options, rho, local_evals, local_iters, stat)
      end if
      result%evaluations = result%evaluations + local_evals
      result%iterations = result%iterations + local_iters
      call evaluate_penalty(problem, candidate, rho, value, grad, violation, .false., stat, result%evaluations)
      if (stat == 0 .and. value < best_value) then
        best_value = value
        best = candidate
      end if
      if (result%evaluations >= options%maxeval) exit
    end do

    call final_values(problem, best, final_f, final_violation, result%evaluations, stat)
    allocate(result%solution(n))
    result%solution = best
    result%objective = final_f
    result%max_constraint = final_violation
    if (result%evaluations >= options%maxeval) then
      result%status = NLOPT_MAXEVAL_REACHED
    else
      result%status = NLOPT_XTOL_REACHED
    end if
    result%converged = final_violation <= max(options%constraint_tol, 1.0e-6_dp)
    result%message = status_message(result%status)
  end subroutine global_search
end module nloptr_solvers
