! SPDX-License-Identifier: GPL-3.0-only
module rgenoud_core
  use rgenoud_kinds, only : dp
  use rgenoud_types, only : genoud_options, genoud_result, objective_fn, &
    gradient_fn, lexical_objective_fn, lexical_better_fn
  use rgenoud_random, only : seed_rng, randu, randi
  use rgenoud_operators, only : oper_uniform, oper_boundary, oper_nonuniform, &
    oper_polytope, oper_simple, oper_whole_nonuniform, oper_heuristic, in_bounds
  use rgenoud_derivatives, only : bfgs_optimize
  implicit none
  private
  public :: run_genoud_core
contains
  subroutine dispatch_lexical_eval(callback, x, f)
    procedure(lexical_objective_fn) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f(:)

    call callback(x, f)
  end subroutine dispatch_lexical_eval

  subroutine run_genoud_core(evalfit, local_fn, local_grad, better, lower, upper, &
      nfit, options, result, local_enabled, gradient_enabled, starting_values)
    procedure(lexical_objective_fn) :: evalfit
    procedure(objective_fn) :: local_fn
    procedure(gradient_fn) :: local_grad
    procedure(lexical_better_fn) :: better
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: nfit
    type(genoud_options), intent(in) :: options
    type(genoud_result), intent(out) :: result
    logical, intent(in) :: local_enabled, gradient_enabled
    real(dp), intent(in), optional :: starting_values(:, :)

    type(genoud_options) :: opt
    real(dp), allocatable :: pop(:, :), newpop(:, :), fit(:, :), newfit(:, :)
    real(dp), allocatable :: cache_x(:, :), cache_f(:, :)
    real(dp), allocatable :: oldbest(:), child(:), p1(:), p2(:), pmat(:, :)
    real(dp), allocatable :: xlocal(:), glocal(:)
    integer, allocatable :: live(:)
    integer :: pop_size, n, generation, max_gen, wait_gen, nochange
    integer :: peak_gen, die_now, op, p2use, i, j, k, cache_n, cache_max
    integer :: targets(9), produced(9), local_iter, eval_requests, unique_evals
    integer :: max_extension_guard
    logical :: accepted, gradient_trigger
    real(dp) :: local_value, mix

    opt = options
    n = size(lower)
    if (size(upper) /= n .or. n < 1 .or. nfit < 1) error stop "invalid genoud dimensions"
    if (any(lower > upper)) error stop "lower bounds exceed upper bounds"
    if (opt%integer_parameters) then
      opt%use_bfgs = .false.
      opt%gradient_check = .false.
      opt%operator_weights(9) = 0.0_dp
    end if

    call seed_rng(opt%seed)
    call set_operator_counts(opt, targets, pop_size)
    max_gen = max(1, opt%max_generations)
    wait_gen = max(0, opt%wait_generations)
    max_extension_guard = max(1000, 16 * max_gen + 16 * wait_gen)

    allocate(pop(n, pop_size), newpop(n, pop_size))
    allocate(fit(nfit, pop_size), newfit(nfit, pop_size))
    allocate(oldbest(nfit), child(n), p1(n), p2(n), live(pop_size))
    allocate(xlocal(n), glocal(n))
    p2use = max(2, n)
    allocate(pmat(n, p2use))

    cache_max = 1
    if (opt%memory_matrix) then
      cache_max = min(opt%max_memory_evals, max(pop_size, (max_gen + 1) * pop_size))
      cache_max = max(cache_max, pop_size)
    end if
    allocate(cache_x(n, cache_max), cache_f(nfit, cache_max))
    cache_n = 0
    eval_requests = 0
    unique_evals = 0

    call initialize_population(pop, lower, upper, opt%integer_parameters, starting_values)
    do j = 1, pop_size
      call evaluate_cached(pop(:, j), fit(:, j))
    end do
    call sort_population(pop, fit, better)
    oldbest = fit(:, 1)
    peak_gen = 0
    nochange = 0
    gradient_trigger = .false.
    generation = 0

    do
      generation = generation + 1
      newpop = pop
      newfit = fit
      produced = 0
      live = 0
      call choose_live(live, sum(targets(2:9)))
      die_now = pop_size

      do while (sum(produced(2:9)) < sum(targets(2:9)))
        op = randi(2, 9)
        if (produced(op) >= targets(op)) cycle

        select case (op)
        case (2)
          i = choose_parent(live)
          if (live(i) > 0) live(i) = live(i) - 1
          child = pop(:, i)
          do k = 1, randi(1, n)
            call oper_uniform(child, lower, upper, opt%integer_parameters)
          end do
          call place_child(child, die_now, newpop, produced(op), 2)

        case (3)
          i = choose_parent(live)
          if (live(i) > 0) live(i) = live(i) - 1
          child = pop(:, i)
          call oper_boundary(child, lower, upper, opt%integer_parameters)
          call place_child(child, die_now, newpop, produced(op), 3)

        case (4)
          i = choose_parent(live)
          if (live(i) > 0) live(i) = live(i) - 1
          child = pop(:, i)
          do k = 1, randi(1, n)
            call oper_nonuniform(child, lower, upper, max_gen, generation, 6, &
              opt%integer_parameters)
          end do
          call place_child(child, die_now, newpop, produced(op), 4)

        case (5)
          do k = 1, p2use
            i = choose_parent(live)
            pmat(:, k) = pop(:, i)
          end do
          call oper_polytope(pmat, child, lower, upper, opt%integer_parameters)
          call place_child(child, die_now, newpop, produced(op), 5)

        case (6)
          if (targets(op) - produced(op) < 2) then
            produced(op) = targets(op)
            cycle
          end if
          call choose_distinct_parents(pop, live, i, j)
          if (live(i) > 0) live(i) = live(i) - 1
          if (live(j) > 0) live(j) = live(j) - 1
          p1 = pop(:, i)
          p2 = pop(:, j)
          call oper_simple(p1, p2, lower, upper, opt%integer_parameters)
          call place_pair(p1, p2, die_now, newpop, produced(op))

        case (7)
          i = choose_parent(live)
          if (live(i) > 0) live(i) = live(i) - 1
          child = pop(:, i)
          call oper_whole_nonuniform(child, lower, upper, max_gen, generation, 6, &
            opt%integer_parameters)
          call place_child(child, die_now, newpop, produced(op), 7)

        case (8)
          if (targets(op) - produced(op) < 2) then
            produced(op) = targets(op)
            cycle
          end if
          call choose_distinct_parents(pop, live, i, j)
          if (live(i) > 0) live(i) = live(i) - 1
          if (live(j) > 0) live(j) = live(j) - 1
          if (better(fit(:, i), fit(:, j))) then
            p1 = pop(:, j)
            p2 = pop(:, i)
          else
            p1 = pop(:, i)
            p2 = pop(:, j)
          end if
          call oper_heuristic(p1, p2, lower, upper, opt%integer_parameters)
          child = p1
          call oper_heuristic(child, p2, lower, upper, opt%integer_parameters)
          call place_pair(p1, child, die_now, newpop, produced(op))

        case (9)
          if (.not. local_enabled .or. generation <= opt%bfgs_burnin) then
            produced(op) = targets(op)
            cycle
          end if
          i = choose_parent(live)
          if (live(i) > 0) live(i) = live(i) - 1
          child = pop(:, i)
          xlocal = child
          call bfgs_optimize(local_fn, local_grad, xlocal, opt%maximize, lower, upper, &
            opt%boundary_enforcement == 2, opt%bfgs_max_iter, opt%bfgs_gtol, &
            local_value, local_iter)
          if (opt%p9_mix > 0.0_dp) then
            mix = min(opt%p9_mix, 1.0_dp)
          else
            mix = randu()
          end if
          accepted = .false.
          do k = 1, 20
            p1 = mix * xlocal + (1.0_dp - mix) * child
            if (opt%boundary_enforcement == 0 .or. in_bounds(p1, lower, upper)) then
              child = p1
              accepted = .true.
              exit
            end if
            mix = 0.5_dp * mix
          end do
          if (.not. accepted) child = pop(:, i)
          call place_child(child, die_now, newpop, produced(op), 9)
        end select
      end do

      pop = newpop
      do j = 1, pop_size
        call evaluate_cached(pop(:, j), fit(:, j))
      end do
      call sort_population(pop, fit, better)

      if (opt%use_bfgs .and. local_enabled .and. generation > opt%bfgs_burnin) then
        xlocal = pop(:, 1)
        call bfgs_optimize(local_fn, local_grad, xlocal, opt%maximize, lower, upper, &
          opt%boundary_enforcement == 2, opt%bfgs_max_iter, opt%bfgs_gtol, &
          local_value, local_iter)
        accepted = opt%boundary_enforcement == 0 .or. in_bounds(xlocal, lower, upper)
        if (accepted) then
          call evalfit(xlocal, newfit(:, 1))
          eval_requests = eval_requests + 1
          unique_evals = unique_evals + 1
          if (better(newfit(:, 1), fit(:, 1))) then
            pop(:, 1) = xlocal
            fit(:, 1) = newfit(:, 1)
          end if
        end if
      end if

      if (significant_improvement(oldbest, fit(:, 1), opt%solution_tolerance, opt%maximize)) then
        oldbest = fit(:, 1)
        nochange = 0
        peak_gen = generation
      else
        nochange = nochange + 1
      end if

      if (opt%print_level > 0) then
        write(*, '(a,i0,a,*(1x,es15.7))') 'generation ', generation, ' best:', fit(:, 1)
      end if

      gradient_trigger = .false.
      if (nochange > wait_gen) then
        if (.not. opt%gradient_check .or. .not. gradient_enabled) then
          exit
        else
          call local_grad(pop(:, 1), glocal)
          if (maxval(abs(glocal)) > opt%solution_tolerance) then
            gradient_trigger = .true.
            wait_gen = max(wait_gen + max(1, wait_gen), wait_gen + 1)
          else
            exit
          end if
        end if
      end if

      if (generation >= max_gen) then
        if (opt%hard_generation_limit) exit
        if (gradient_trigger .or. nochange < wait_gen) then
          max_gen = max_gen + max(max_gen, wait_gen)
        else
          exit
        end if
      end if
      if (generation >= max_extension_guard) exit
    end do

    allocate(result%par(n), result%fit(nfit), result%gradient(n))
    result%par = pop(:, 1)
    result%fit = fit(:, 1)
    result%gradient = 0.0_dp
    if (gradient_enabled) call local_grad(result%par, result%gradient)
    result%generations = generation
    result%peak_generation = peak_gen
    result%pop_size = pop_size
    result%operators = targets
    result%evaluations = eval_requests
    result%unique_evaluations = unique_evals
    result%converged = nochange > wait_gen .and. .not. gradient_trigger
    result%status = merge(0, 1, result%converged)

  contains
    subroutine evaluate_cached(x, fval)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: fval(:)
      integer :: q
      eval_requests = eval_requests + 1
      if (opt%memory_matrix) then
        do q = 1, cache_n
          if (all(abs(x - cache_x(:, q)) <= 0.0_dp)) then
            fval = cache_f(:, q)
            return
          end if
        end do
      end if
      call dispatch_lexical_eval(evalfit, x, fval)
      unique_evals = unique_evals + 1
      if (opt%memory_matrix .and. cache_n < cache_max) then
        cache_n = cache_n + 1
        cache_x(:, cache_n) = x
        cache_f(:, cache_n) = fval
      end if
    end subroutine evaluate_cached

    subroutine choose_live(counts, needed)
      integer, intent(out) :: counts(:)
      integer, intent(in) :: needed
      integer :: q, idx
      counts = 0
      do q = 1, needed
        idx = rank_draw(size(counts))
        counts(idx) = counts(idx) + 1
      end do
    end subroutine choose_live

    integer function rank_draw(m) result(idx)
      integer, intent(in) :: m
      real(dp) :: u, c, p
      integer :: q
      u = randu()
      c = 0.0_dp
      p = 0.5_dp
      idx = m
      do q = 1, m
        c = c + p
        if (u <= c .or. q == m) then
          idx = q
          return
        end if
        p = 0.5_dp * p
      end do
    end function rank_draw

    integer function choose_parent(counts) result(idx)
      integer, intent(in) :: counts(:)
      integer :: total, pick, q, acc
      total = sum(max(counts, 0))
      if (total <= 0) then
        idx = rank_draw(size(counts))
        return
      end if
      pick = randi(1, total)
      acc = 0
      do q = 1, size(counts)
        acc = acc + max(counts(q), 0)
        if (pick <= acc) then
          idx = q
          return
        end if
      end do
      idx = 1
    end function choose_parent

    subroutine choose_distinct_parents(xp, counts, ia, ib)
      real(dp), intent(in) :: xp(:, :)
      integer, intent(in) :: counts(:)
      integer, intent(out) :: ia, ib
      integer :: tries
      ia = choose_parent(counts)
      ib = choose_parent(counts)
      do tries = 1, 1000
        if (ia /= ib .and. any(abs(xp(:, ia) - xp(:, ib)) > 0.0_dp)) return
        ib = choose_parent(counts)
      end do
      ib = min(size(xp, 2), ia + 1)
      if (ib == ia) ib = max(1, ia - 1)
    end subroutine choose_distinct_parents

    subroutine place_child(x, dnow, xnew, made, operator_number)
      real(dp), intent(in) :: x(:)
      integer, intent(inout) :: dnow, made
      real(dp), intent(inout) :: xnew(:, :)
      integer, intent(in) :: operator_number
      if (dnow < 2) then
        made = targets(operator_number)
        return
      end if
      xnew(:, dnow) = x
      dnow = dnow - 1
      made = made + 1
    end subroutine place_child

    subroutine place_pair(a, b, dnow, xnew, made)
      real(dp), intent(in) :: a(:), b(:)
      integer, intent(inout) :: dnow, made
      real(dp), intent(inout) :: xnew(:, :)
      if (dnow < 3) then
        made = made + 2
        return
      end if
      xnew(:, dnow) = a
      dnow = dnow - 1
      xnew(:, dnow) = b
      dnow = dnow - 1
      made = made + 2
    end subroutine place_pair
  end subroutine run_genoud_core

  subroutine set_operator_counts(opt, counts, pop_size)
    type(genoud_options), intent(in) :: opt
    integer, intent(out) :: counts(9), pop_size
    real(dp) :: weights(9), total
    integer :: ptotal, i

    weights = max(opt%operator_weights, 0.0_dp)
    if (opt%integer_parameters) weights(9) = 0.0_dp
    pop_size = max(4, opt%pop_size)
    total = sum(weights)
    counts = 0
    if (total > 0.0_dp) then
      do i = 1, 9
        counts(i) = nint((weights(i) / total) * real(pop_size - 2, dp))
      end do
    end if
    if (mod(counts(6), 2) /= 0) counts(6) = counts(6) + 1
    if (mod(counts(8), 2) /= 0) counts(8) = counts(8) + 1
    ptotal = sum(counts(2:9))
    if (ptotal >= pop_size) then
      pop_size = ptotal + 1
      if (mod(pop_size, 2) /= 0) pop_size = pop_size + 1
    else if (mod(pop_size, 2) /= 0) then
      pop_size = pop_size + 1
    end if
    counts(1) = pop_size - ptotal - 1
  end subroutine set_operator_counts

  subroutine initialize_population(pop, lower, upper, integer_mode, starting_values)
    real(dp), intent(out) :: pop(:, :)
    real(dp), intent(in) :: lower(:), upper(:)
    logical, intent(in) :: integer_mode
    real(dp), intent(in), optional :: starting_values(:, :)
    integer :: i, j, ns, lo, hi

    ns = 0
    if (present(starting_values)) then
      if (size(starting_values, 1) /= size(pop, 1)) then
        error stop "starting_values must have nvars rows"
      end if
      ns = min(size(starting_values, 2), size(pop, 2))
      pop(:, 1:ns) = starting_values(:, 1:ns)
      if (integer_mode) pop(:, 1:ns) = real(nint(pop(:, 1:ns)), dp)
    end if
    do j = ns + 1, size(pop, 2)
      do i = 1, size(pop, 1)
        if (integer_mode) then
          lo = ceiling(lower(i))
          hi = floor(upper(i))
          pop(i, j) = real(randi(lo, hi), dp)
        else
          pop(i, j) = randu(lower(i), upper(i))
        end if
      end do
    end do
  end subroutine initialize_population

  subroutine sort_population(pop, fit, better)
    real(dp), intent(inout) :: pop(:, :), fit(:, :)
    procedure(lexical_better_fn) :: better
    real(dp) :: xtmp(size(pop, 1)), ftmp(size(fit, 1))
    integer :: i, j

    do i = 2, size(pop, 2)
      xtmp = pop(:, i)
      ftmp = fit(:, i)
      j = i - 1
      do while (j >= 1)
        if (.not. better(ftmp, fit(:, j))) exit
        pop(:, j + 1) = pop(:, j)
        fit(:, j + 1) = fit(:, j)
        j = j - 1
      end do
      pop(:, j + 1) = xtmp
      fit(:, j + 1) = ftmp
    end do
  end subroutine sort_population

  logical function significant_improvement(oldf, newf, tol, maximize) result(ok)
    real(dp), intent(in) :: oldf(:), newf(:), tol
    logical, intent(in) :: maximize
    if (maximize) then
      ok = any(newf > oldf + tol)
    else
      ok = any(newf < oldf - tol)
    end if
  end function significant_improvement
end module rgenoud_core
