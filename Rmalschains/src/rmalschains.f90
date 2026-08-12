module rmalschains
  use iso_fortran_env, only : int64
  use rmalschains_kinds, only : dp
  use rmalschains_rng, only : rng_state, rng_seed, rng_uniform, rng_int
  use rmalschains_types, only : mals_control, mals_result, mals_objective, ls_state
  use rmalschains_operators, only : blx_alpha, bga_mutate, nam_select
  use rmalschains_local_search, only : apply_local_search, reset_ls_state
  implicit none
  private
  public :: dp, mals_control, mals_result, mals_objective, malschains_optimize, malschains_control
contains
  function malschains_control(popsize, ls, istep, effort, alpha, optimum, threshold, ls_only, &
    ls_param1, ls_param2, seed) result(ctrl)
    integer, intent(in), optional :: popsize, istep, ls_param1, ls_param2
    character(len=*), intent(in), optional :: ls
    real(dp), intent(in), optional :: effort, alpha, optimum, threshold
    logical, intent(in), optional :: ls_only
    integer(int64), intent(in), optional :: seed
    type(mals_control) :: ctrl
    integer :: p
    ctrl = mals_control()
    if (present(popsize)) then
      p = nint(real(popsize, dp) / 10.0_dp) * 10
      if (p == 0) p = 10
      ctrl%popsize = p
    end if
    if (present(ls)) ctrl%ls = trim(ls)
    if (present(istep)) ctrl%istep = nint(real(istep, dp))
    if (present(effort)) ctrl%effort = effort
    if (present(alpha)) ctrl%alpha = alpha
    if (present(optimum)) ctrl%optimum = optimum
    if (present(threshold)) ctrl%threshold = threshold
    if (present(ls_only)) ctrl%ls_only = ls_only
    if (present(ls_param1)) ctrl%ls_param1 = ls_param1
    if (present(ls_param2)) ctrl%ls_param2 = ls_param2
    if (present(seed)) ctrl%seed = seed
  end function malschains_control

  function malschains_optimize(fn, lower, upper, max_evals, control, initialpop) result(res)
    procedure(mals_objective) :: fn
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: max_evals
    type(mals_control), intent(in), optional :: control
    real(dp), intent(in), optional :: initialpop(:, :)
    type(mals_result) :: res

    type(mals_control) :: ctrl
    type(rng_state) :: rng
    type(ls_state), allocatable :: states(:)
    real(dp), allocatable :: pop(:, :), fit(:), child(:), x(:)
    integer, allocatable :: non_improved(:)
    integer :: n, psize, i, j, ninit, mom, dad, worst, best, pos, ea_budget, used
    integer :: remaining, old_best_idx, new_best_idx
    real(dp) :: fx, oldfx, old_best, new_best, old_second, new_second, fchild
    real(dp) :: t_ma0, t0, t1
    logical :: ea_improved, ls_improved, mut_changed

    ctrl = mals_control()
    if (present(control)) ctrl = control
    call validate_control(ctrl)
    n = size(lower)
    if (n < 1 .or. size(upper) /= n) error stop "malschains_optimize: invalid bounds"
    if (any(lower >= upper)) error stop "malschains_optimize: lower must be below upper"
    if (max_evals < 1) error stop "malschains_optimize: max_evals must be positive"
    psize = nint(real(ctrl%popsize, dp) / 10.0_dp) * 10
    if (psize == 0) psize = 10
    if (psize < 10) psize = 10
    call rng_seed(rng, ctrl%seed)
    allocate(pop(n, psize), fit(psize), non_improved(psize), states(psize), child(n), x(n))
    non_improved = 0
    fit = huge(1.0_dp)
    do j = 1, psize
      do i = 1, n
        pop(i, j) = lower(i) + rng_uniform(rng) * (upper(i) - lower(i))
      end do
    end do
    call cpu_time(t_ma0)

    if (ctrl%ls_only) then
      call run_ls_only(fn, lower, upper, max_evals, ctrl, initialpop, pop, fit, states, rng, res)
      call cpu_time(t1)
      res%time_ms_ma = 1000.0_dp * (t1 - t_ma0)
      return
    end if

    call cpu_time(t0)
    do j = 1, psize
      fit(j) = fn(pop(:, j))
      res%actual_nfe = res%actual_nfe + 1
      res%num_eval_ea = res%num_eval_ea + 1
    end do
    call cpu_time(t1)
    res%time_ms_ea = res%time_ms_ea + 1000.0_dp * (t1 - t0)

    if (present(initialpop)) then
      if (size(initialpop, 1) /= n) error stop "malschains_optimize: initialpop rows must equal dimension"
      ninit = min(size(initialpop, 2), psize)
      do j = 1, ninit
        pop(:, j) = initialpop(:, j)
        if (any(pop(:, j) < lower) .or. any(pop(:, j) > upper)) &
          error stop "malschains_optimize: initial population outside bounds"
        fit(j) = fn(pop(:, j))
        res%actual_nfe = res%actual_nfe + 1
      end do
    end if

    do while (res%num_eval_ea + res%num_eval_ls < max_evals)
      best = minloc(fit, dim=1)
      if (target_hit(fit(best), ctrl)) exit
      old_best_idx = best
      old_best = fit(best)
      old_second = second_best(fit)
      remaining = max_evals - (res%num_eval_ea + res%num_eval_ls)
      ea_budget = calculate_frec(res%num_eval_ea, res%num_eval_ls, ctrl%istep, ctrl%effort)
      ea_budget = min(max(ea_budget, 0), remaining)

      call cpu_time(t0)
      do i = 1, ea_budget
        if (target_hit(minval(fit), ctrl)) exit
        call nam_select(pop, rng, mom, dad)
        call blx_alpha(pop(:, mom), pop(:, dad), lower, upper, ctrl%alpha, rng, child)
        call bga_mutate(child, lower, upper, rng, mut_changed)
        fchild = fn(child)
        res%actual_nfe = res%actual_nfe + 1
        res%num_eval_ea = res%num_eval_ea + 1
        worst = maxloc(fit, dim=1)
        if (fchild < fit(worst)) then
          pop(:, worst) = child
          fit(worst) = fchild
          non_improved(worst) = 0
          call reset_ls_state(states(worst))
        end if
      end do
      call cpu_time(t1)
      res%time_ms_ea = res%time_ms_ea + 1000.0_dp * (t1 - t0)
      res%generations_ea = res%generations_ea + ea_budget
      new_best_idx = minloc(fit, dim=1)
      new_best = fit(new_best_idx)
      new_second = second_best(fit)
      res%improvement_ea = res%improvement_ea + abs(new_best - old_best)
      if (new_best < old_best) res%num_improvement_ea = res%num_improvement_ea + 1
      res%num_total_ea = res%num_total_ea + 1
      ea_improved = improved_enough(old_second, new_second, ctrl%threshold)
      if (target_hit(new_best, ctrl)) exit
      remaining = max_evals - (res%num_eval_ea + res%num_eval_ls)
      if (remaining <= 0) exit

      pos = choose_ls_individual(fit, non_improved, rng)
      x = pop(:, pos)
      fx = fit(pos)
      oldfx = fx
      call cpu_time(t0)
      call apply_local_search(ctrl%ls, states(pos), x, fx, min(ctrl%istep, remaining), pop, pos, &
        lower, upper, rng, fn, ctrl, used)
      call cpu_time(t1)
      res%time_ms_ls = res%time_ms_ls + 1000.0_dp * (t1 - t0)
      res%actual_nfe = res%actual_nfe + used
      res%num_eval_ls = res%num_eval_ls + used
      res%local_search_calls = res%local_search_calls + 1
      res%improvement_ls = res%improvement_ls + abs(fx - oldfx)
      ls_improved = improved_enough(new_best, fx, ctrl%threshold)
      if (fx < new_best) res%num_improvement_ls = res%num_improvement_ls + 1
      res%num_total_ls = res%num_total_ls + 1
      if (improved_enough(oldfx, fx, ctrl%threshold)) then
        pop(:, pos) = x
        fit(pos) = fx
      else
        non_improved(pos) = non_improved(pos) + 1
        call reset_ls_state(states(pos))
      end if
      if (.not. ea_improved .and. .not. ls_improved) then
        continue
      end if
    end do

    best = minloc(fit, dim=1)
    allocate(res%sol(n))
    res%sol = pop(:, best)
    res%fitness = fit(best)
    call cpu_time(t1)
    res%time_ms_ma = 1000.0_dp * (t1 - t_ma0)
    if (target_hit(res%fitness, ctrl)) then
      res%stop_reason = 'target reached'
    else if (res%num_eval_ea + res%num_eval_ls >= max_evals) then
      res%stop_reason = 'maximum evaluations reached'
    else
      res%stop_reason = 'terminated'
    end if
  end function malschains_optimize

  subroutine run_ls_only(fn, lower, upper, max_evals, ctrl, initialpop, pop, fit, states, rng, res)
    procedure(mals_objective) :: fn
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(in) :: max_evals
    type(mals_control), intent(in) :: ctrl
    real(dp), intent(in), optional :: initialpop(:, :)
    real(dp), intent(inout) :: pop(:, :), fit(:)
    type(ls_state), intent(inout) :: states(:)
    type(rng_state), intent(inout) :: rng
    type(mals_result), intent(inout) :: res
    real(dp), allocatable :: x(:)
    real(dp) :: fx, oldfx, t0, t1
    integer :: j, ninit, pos, used, accounted
    allocate(x(size(lower)))
    if (present(initialpop)) then
      if (size(initialpop, 1) /= size(lower)) error stop "run_ls_only: initialpop dimension mismatch"
      ninit = min(size(initialpop, 2), size(pop, 2))
      if (ninit < 1) error stop "run_ls_only: initialpop is empty"
      do j = 1, ninit
        pop(:, j) = initialpop(:, j)
        fit(j) = fn(pop(:, j))
        res%actual_nfe = res%actual_nfe + 1
      end do
      pos = minloc(fit(1:ninit), dim=1)
      x = pop(:, pos); fx = fit(pos)
    else if (ctrl%legacy_ls_only_zero_start) then
      x = 0.0_dp
      fx = 0.0_dp
      pos = 1
    else
      do j = 1, size(pop, 2)
        fit(j) = fn(pop(:, j))
        res%actual_nfe = res%actual_nfe + 1
      end do
      pos = minloc(fit, dim=1)
      x = pop(:, pos); fx = fit(pos)
    end if
    accounted = 0
    oldfx = 0.0_dp
    do while (accounted < max_evals .and. .not. target_hit(fx, ctrl))
      call cpu_time(t0)
      call apply_local_search(ctrl%ls, states(pos), x, fx, min(100, max_evals - accounted), &
        pop, pos, lower, upper, rng, fn, ctrl, used)
      call cpu_time(t1)
      res%time_ms_ls = res%time_ms_ls + 1000.0_dp * (t1 - t0)
      res%actual_nfe = res%actual_nfe + used
      res%local_search_calls = res%local_search_calls + 1
      if (abs(fx - oldfx) < ctrl%threshold) exit
      accounted = accounted + 100
      oldfx = fx
      if (used <= 0) exit
    end do
    res%num_eval_ls = min(accounted, max_evals)
    allocate(res%sol(size(x)))
    res%sol = x; res%fitness = fx
    if (target_hit(fx, ctrl)) then
      res%stop_reason = 'target reached'
    else if (accounted >= max_evals) then
      res%stop_reason = 'maximum evaluations reached'
    else
      res%stop_reason = 'local-search improvement below threshold'
    end if
  end subroutine run_ls_only

  subroutine validate_control(ctrl)
    type(mals_control), intent(in) :: ctrl
    if (ctrl%popsize <= 0) error stop "malschains: popsize must be positive"
    if (ctrl%istep <= 0) error stop "malschains: istep must be positive"
    if (.not. ctrl%ls_only) then
      if (ctrl%effort <= 0.0_dp .or. ctrl%effort >= 1.0_dp) &
        error stop "malschains: effort must lie strictly between zero and one"
    end if
    if (ctrl%alpha <= 0.0_dp) error stop "malschains: alpha must be positive"
    if (ctrl%threshold < 0.0_dp) error stop "malschains: threshold must be nonnegative"
  end subroutine validate_control

  pure function target_hit(fx, ctrl) result(hit)
    real(dp), intent(in) :: fx
    type(mals_control), intent(in) :: ctrl
    logical :: hit
    hit = fx <= ctrl%optimum + ctrl%threshold
  end function target_hit

  pure function improved_enough(oldfx, newfx, threshold) result(ok)
    real(dp), intent(in) :: oldfx, newfx, threshold
    logical :: ok
    real(dp) :: minimum
    minimum = threshold / 10.0_dp
    if (newfx >= oldfx) then
      ok = .false.
    else if (abs(newfx - oldfx) < minimum) then
      ok = .false.
    else if (threshold <= 0.0_dp) then
      ok = .true.
    else
      ok = abs(oldfx - newfx) / max(abs(newfx), tiny(1.0_dp)) >= threshold
    end if
  end function improved_enough

  pure function calculate_frec(neval_ea, neval_ls, intensity, ratio) result(frec)
    integer, intent(in) :: neval_ea, neval_ls, intensity
    real(dp), intent(in) :: ratio
    integer :: frec
    real(dp) :: value
    value = (ratio * real(neval_ea + neval_ls + intensity, dp) - real(neval_ea, dp)) / (1.0_dp - ratio)
    frec = int(floor(value))
  end function calculate_frec

  function choose_ls_individual(fit, non_improved, rng) result(pos)
    real(dp), intent(in) :: fit(:)
    integer, intent(in) :: non_improved(:)
    type(rng_state), intent(inout) :: rng
    integer :: pos, i
    real(dp) :: best
    pos = 0; best = huge(1.0_dp)
    do i = 1, size(fit)
      if (non_improved(i) == 0 .and. fit(i) < best) then
        best = fit(i); pos = i
      end if
    end do
    if (pos == 0) pos = rng_int(rng, 1, size(fit))
  end function choose_ls_individual

  pure function second_best(fit) result(value)
    real(dp), intent(in) :: fit(:)
    real(dp) :: value, best
    integer :: i, ibest
    ibest = minloc(fit, dim=1)
    best = fit(ibest)
    value = huge(1.0_dp)
    do i = 1, size(fit)
      if (i /= ibest .and. fit(i) < value) value = fit(i)
    end do
    if (size(fit) == 1) value = best
  end function second_best
end module rmalschains
