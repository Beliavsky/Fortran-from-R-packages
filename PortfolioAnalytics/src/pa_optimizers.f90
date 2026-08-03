! SPDX-License-Identifier: GPL-3.0-only
module pa_optimizers
  use pa_kinds, only: dp, pa_huge
  use pa_types, only: portfolio_constraints, portfolio_options, portfolio_result, &
       opt_auto, opt_projected_gradient, opt_differential_evolution, opt_random_search, &
       obj_min_variance, obj_max_return, obj_quadratic_utility, pa_success, &
       pa_invalid_input, pa_infeasible, pa_max_iterations
  use pa_linalg, only: set_random_seed
  use pa_statistics, only: sample_moments, portfolio_mean, portfolio_stddev, &
       portfolio_returns, historical_es, turnover, hhi, diversification, risk_contributions
  use pa_constraints, only: repair_weights, check_constraints, random_feasible_portfolio
  use pa_objectives, only: evaluate_portfolio_objective, objective_gradient_numeric
  implicit none
  private
  public :: optimize_portfolio, equal_weight_portfolio, inverse_volatility_portfolio

contains

  subroutine equal_weight_portfolio(c, weights, feasible)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(out) :: weights(:)
    logical, intent(out) :: feasible
    real(dp), allocatable :: candidate(:)
    integer :: n
    n = size(weights)
    allocate(candidate(n))
    candidate = 0.5_dp*(c%min_sum+c%max_sum)/real(n,dp)
    call repair_weights(candidate,c,weights,feasible)
  end subroutine equal_weight_portfolio

  subroutine inverse_volatility_portfolio(sigma, c, weights, feasible)
    real(dp), intent(in) :: sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(out) :: weights(:)
    logical, intent(out) :: feasible
    real(dp), allocatable :: candidate(:)
    real(dp) :: target
    integer :: i, n
    n = size(weights)
    allocate(candidate(n))
    do i = 1, n
      candidate(i) = 1.0_dp/sqrt(max(sigma(i,i),1.0e-16_dp))
    end do
    target = 0.5_dp*(c%min_sum+c%max_sum)
    candidate = target*candidate/sum(candidate)
    call repair_weights(candidate,c,weights,feasible)
  end subroutine inverse_volatility_portfolio

  subroutine optimize_portfolio(returns, c, options, result, mu_input, sigma_input, budgets)
    real(dp), intent(in) :: returns(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    type(portfolio_result), intent(out) :: result
    real(dp), intent(in), optional :: mu_input(:), sigma_input(:,:), budgets(:)
    real(dp), allocatable :: mu(:), sigma(:,:), weights(:)
    integer :: n, method, evals, iters
    logical :: converged

    n = size(returns,2)
    result%status = pa_invalid_input
    result%message = 'invalid input'
    if (size(returns,1) < 2 .or. n < 1) return
    if (.not. allocated(c%min_weight) .or. .not. allocated(c%max_weight)) return
    if (size(c%min_weight) /= n .or. size(c%max_weight) /= n) return
    allocate(mu(n),sigma(n,n),weights(n))
    if (present(mu_input)) then
      if (size(mu_input) /= n) return
      mu = mu_input
    else
      call sample_moments(returns,mu,sigma)
    end if
    if (present(sigma_input)) then
      if (size(sigma_input,1) /= n .or. size(sigma_input,2) /= n) return
      sigma = sigma_input
    else if (present(mu_input)) then
      call sample_moments(returns,weights,sigma)
    end if

    method = options%optimizer
    if (method == opt_auto) then
      select case(options%objective)
      case(obj_min_variance,obj_max_return,obj_quadratic_utility)
        method = opt_projected_gradient
      case default
        method = opt_differential_evolution
      end select
    end if

    select case(method)
    case(opt_projected_gradient)
      if (present(budgets)) then
        call projected_gradient_solver(returns,mu,sigma,c,options,weights,converged,iters,evals,budgets)
      else
        call projected_gradient_solver(returns,mu,sigma,c,options,weights,converged,iters,evals)
      end if
    case(opt_differential_evolution)
      if (present(budgets)) then
        call differential_evolution_solver(returns,mu,sigma,c,options,weights,converged,iters,evals,budgets)
      else
        call differential_evolution_solver(returns,mu,sigma,c,options,weights,converged,iters,evals)
      end if
    case(opt_random_search)
      if (present(budgets)) then
        call random_search_solver(returns,mu,sigma,c,options,weights,converged,iters,evals,budgets)
      else
        call random_search_solver(returns,mu,sigma,c,options,weights,converged,iters,evals)
      end if
    case default
      return
    end select

    if (options%local_refine .and. method /= opt_projected_gradient) then
      if (present(budgets)) then
        call pattern_refine(returns,mu,sigma,c,options,weights,evals,budgets)
      else
        call pattern_refine(returns,mu,sigma,c,options,weights,evals)
      end if
    end if
    if (present(budgets)) then
      call fill_result(returns,mu,sigma,c,options,weights,converged,iters,evals,result,budgets)
    else
      call fill_result(returns,mu,sigma,c,options,weights,converged,iters,evals,result)
    end if
  end subroutine optimize_portfolio

  subroutine initial_weights(c, mu, seed, weights)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(in) :: mu(:)
    integer, intent(in) :: seed
    real(dp), intent(out) :: weights(:)
    real(dp), allocatable :: candidate(:)
    logical :: ok
    integer :: n
    n = size(weights)
    allocate(candidate(n))
    if (allocated(c%initial_weights)) then
      candidate = c%initial_weights
    else
      candidate = 0.5_dp*(c%min_sum+c%max_sum)/real(n,dp)
    end if
    call repair_weights(candidate,c,weights,ok)
    if (.not. check_constraints(weights,c,mu)) then
      call random_feasible_portfolio(c,weights,ok,seed,mu,20000)
      if (.not. ok) then
        candidate = 0.5_dp*(c%min_weight+c%max_weight)
        call repair_weights(candidate,c,weights,ok)
      end if
    end if
  end subroutine initial_weights

  subroutine projected_gradient_solver(returns,mu,sigma,c,options,best,converged,iters,evals,budgets)
    real(dp), intent(in) :: returns(:,:),mu(:),sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(out) :: best(:)
    logical, intent(out) :: converged
    integer, intent(out) :: iters,evals
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: grad(:), candidate(:), trial(:)
    real(dp) :: f, ftrial, step, gnorm, change
    logical :: ok
    integer :: iter, ls, n

    n = size(best)
    allocate(grad(n),candidate(n),trial(n))
    call initial_weights(c,mu,options%seed,best)
    if (present(budgets)) then
      f = evaluate_portfolio_objective(best,returns,mu,sigma,c,options,budgets)
    else
      f = evaluate_portfolio_objective(best,returns,mu,sigma,c,options)
    end if
    evals = 1
    converged = .false.
    do iter = 1, options%max_iterations
      select case(options%objective)
      case(obj_min_variance)
        grad = 2.0_dp*matmul(sigma,best)
        evals = evals + 1
      case(obj_max_return)
        grad = -mu
        evals = evals + 1
      case(obj_quadratic_utility)
        grad = options%risk_aversion*matmul(sigma,best)-mu
        evals = evals + 1
      case default
        if (present(budgets)) then
          call objective_gradient_numeric(best,returns,mu,sigma,c,options,grad,budgets)
        else
          call objective_gradient_numeric(best,returns,mu,sigma,c,options,grad)
        end if
        evals = evals + 2*n
      end select
      if (abs(options%concentration_aversion)>0.0_dp) &
        grad=grad+2.0_dp*options%concentration_aversion*best
      if (allocated(c%initial_weights) .and. abs(options%turnover_aversion)>0.0_dp) &
        grad=grad+options%turnover_aversion*sign(1.0_dp,best-c%initial_weights)/real(n,dp)
      gnorm = sqrt(dot_product(grad,grad))
      if (gnorm <= options%tolerance) then
        converged = .true.
        exit
      end if
      step = 1.0_dp/max(gnorm,1.0_dp)
      do ls = 1, 30
        candidate = best-step*grad
        call repair_weights(candidate,c,trial,ok)
        if (.not. ok) then
          step = 0.5_dp*step
          cycle
        end if
        if (present(budgets)) then
          ftrial = evaluate_portfolio_objective(trial,returns,mu,sigma,c,options,budgets)
        else
          ftrial = evaluate_portfolio_objective(trial,returns,mu,sigma,c,options)
        end if
        evals = evals + 1
        if (ftrial < f-1.0e-4_dp*step*gnorm*gnorm .or. ftrial < f) exit
        step = 0.5_dp*step
      end do
      if (ls > 30) then
        converged = .true.
        exit
      end if
      change = maxval(abs(trial-best))
      best = trial
      f = ftrial
      if (change <= options%tolerance) then
        converged = .true.
        exit
      end if
    end do
    iters = min(iter,options%max_iterations)
  end subroutine projected_gradient_solver

  subroutine differential_evolution_solver(returns,mu,sigma,c,options,best,converged,iters,evals,budgets)
    real(dp), intent(in) :: returns(:,:),mu(:),sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(out) :: best(:)
    logical, intent(out) :: converged
    integer, intent(out) :: iters,evals
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: pop(:,:), values(:), mutant(:), trial(:), repaired(:)
    real(dp) :: fscale, cr, u, ftrial, best_old
    integer :: n, npop, i, j, r1, r2, r3, jrand, iter, best_idx
    logical :: ok

    n = size(best)
    npop = max(options%population_size,4*n+4)
    allocate(pop(n,npop),values(npop),mutant(n),trial(n),repaired(n))
    call set_random_seed(options%seed)
    call initial_weights(c,mu,options%seed,pop(:,1))
    do i = 2, npop
      call random_feasible_portfolio(c,pop(:,i),ok,options%seed+i,mu,5000)
      if (.not. ok) then
        call random_number(pop(:,i))
        pop(:,i) = c%min_weight + pop(:,i)*(c%max_weight-c%min_weight)
        call repair_weights(pop(:,i),c,repaired,ok)
        if (ok) pop(:,i) = repaired
      end if
    end do
    do i = 1, npop
      if (present(budgets)) then
        values(i) = evaluate_portfolio_objective(pop(:,i),returns,mu,sigma,c,options,budgets)
      else
        values(i) = evaluate_portfolio_objective(pop(:,i),returns,mu,sigma,c,options)
      end if
    end do
    evals = npop
    fscale = 0.75_dp
    cr = 0.90_dp
    converged = .false.
    do iter = 1, options%max_iterations
      best_old = minval(values)
      do i = 1, npop
        call draw_distinct_indices(npop,i,r1,r2,r3)
        mutant = pop(:,r1)+fscale*(pop(:,r2)-pop(:,r3))
        trial = pop(:,i)
        call random_number(u)
        jrand = 1+int(u*real(n,dp))
        jrand = min(max(jrand,1),n)
        do j = 1, n
          call random_number(u)
          if (u <= cr .or. j == jrand) trial(j) = mutant(j)
        end do
        call repair_weights(trial,c,repaired,ok)
        if (.not. ok) cycle
        if (present(budgets)) then
          ftrial = evaluate_portfolio_objective(repaired,returns,mu,sigma,c,options,budgets)
        else
          ftrial = evaluate_portfolio_objective(repaired,returns,mu,sigma,c,options)
        end if
        evals = evals+1
        if (ftrial <= values(i)) then
          pop(:,i) = repaired
          values(i) = ftrial
        end if
      end do
      if (abs(best_old-minval(values)) <= options%tolerance*max(1.0_dp,abs(best_old))) then
        if (population_spread(pop) <= sqrt(options%tolerance)) then
          converged = .true.
          exit
        end if
      end if
    end do
    best_idx = minloc(values,dim=1)
    best = pop(:,best_idx)
    iters = min(iter,options%max_iterations)
  end subroutine differential_evolution_solver

  subroutine draw_distinct_indices(n, excluded, r1, r2, r3)
    integer, intent(in) :: n, excluded
    integer, intent(out) :: r1, r2, r3
    real(dp) :: u
    do
      call random_number(u)
      r1 = min(n,1+int(u*real(n,dp)))
      if (r1 /= excluded) exit
    end do
    do
      call random_number(u)
      r2 = min(n,1+int(u*real(n,dp)))
      if (r2 /= excluded .and. r2 /= r1) exit
    end do
    do
      call random_number(u)
      r3 = min(n,1+int(u*real(n,dp)))
      if (r3 /= excluded .and. r3 /= r1 .and. r3 /= r2) exit
    end do
  end subroutine draw_distinct_indices

  real(dp) function population_spread(pop) result(value)
    real(dp), intent(in) :: pop(:,:)
    real(dp), allocatable :: meanv(:)
    integer :: i
    allocate(meanv(size(pop,1)))
    meanv = sum(pop,dim=2)/real(size(pop,2),dp)
    value = 0.0_dp
    do i = 1, size(pop,2)
      value = max(value,maxval(abs(pop(:,i)-meanv)))
    end do
  end function population_spread

  subroutine random_search_solver(returns,mu,sigma,c,options,best,converged,iters,evals,budgets)
    real(dp), intent(in) :: returns(:,:),mu(:),sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(out) :: best(:)
    logical, intent(out) :: converged
    integer, intent(out) :: iters,evals
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: candidate(:)
    real(dp) :: fbest, f
    integer :: i
    logical :: ok
    allocate(candidate(size(best)))
    call initial_weights(c,mu,options%seed,best)
    if (present(budgets)) then
      fbest = evaluate_portfolio_objective(best,returns,mu,sigma,c,options,budgets)
    else
      fbest = evaluate_portfolio_objective(best,returns,mu,sigma,c,options)
    end if
    evals = 1
    do i = 1, options%random_portfolios
      call random_feasible_portfolio(c,candidate,ok,options%seed+i,mu,5000)
      if (.not. ok) cycle
      if (present(budgets)) then
        f = evaluate_portfolio_objective(candidate,returns,mu,sigma,c,options,budgets)
      else
        f = evaluate_portfolio_objective(candidate,returns,mu,sigma,c,options)
      end if
      evals = evals+1
      if (f < fbest) then
        fbest = f
        best = candidate
      end if
    end do
    iters = options%random_portfolios
    converged = check_constraints(best,c,mu)
  end subroutine random_search_solver

  subroutine pattern_refine(returns,mu,sigma,c,options,weights,evals,budgets)
    real(dp), intent(in) :: returns(:,:),mu(:),sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(inout) :: weights(:)
    integer, intent(inout) :: evals
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: candidate(:), trial(:)
    real(dp) :: f, ft, step
    integer :: i, direction, sweep
    logical :: ok, improved
    allocate(candidate(size(weights)),trial(size(weights)))
    if (present(budgets)) then
      f = evaluate_portfolio_objective(weights,returns,mu,sigma,c,options,budgets)
    else
      f = evaluate_portfolio_objective(weights,returns,mu,sigma,c,options)
    end if
    step = 0.1_dp*max(1.0_dp,maxval(c%max_weight-c%min_weight))
    do sweep = 1, 100
      improved = .false.
      do i = 1, size(weights)
        do direction = -1, 1, 2
          candidate = weights
          candidate(i) = candidate(i)+real(direction,dp)*step
          call repair_weights(candidate,c,trial,ok)
          if (.not. ok) cycle
          if (present(budgets)) then
            ft = evaluate_portfolio_objective(trial,returns,mu,sigma,c,options,budgets)
          else
            ft = evaluate_portfolio_objective(trial,returns,mu,sigma,c,options)
          end if
          evals = evals+1
          if (ft < f) then
            weights = trial
            f = ft
            improved = .true.
          end if
        end do
      end do
      if (.not. improved) step = 0.5_dp*step
      if (step <= options%tolerance) exit
    end do
  end subroutine pattern_refine

  subroutine fill_result(returns,mu,sigma,c,options,weights,converged,iters,evals,result,budgets)
    real(dp), intent(in) :: returns(:,:),mu(:),sigma(:,:),weights(:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    logical, intent(in) :: converged
    integer, intent(in) :: iters,evals
    type(portfolio_result), intent(out) :: result
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: rp(:)
    integer :: n
    n = size(weights)
    allocate(result%weights(n),result%risk_contribution(n),rp(size(returns,1)))
    result%weights = weights
    if (present(budgets)) then
      result%objective_value = evaluate_portfolio_objective(weights,returns,mu,sigma,c,options,budgets)
    else
      result%objective_value = evaluate_portfolio_objective(weights,returns,mu,sigma,c,options)
    end if
    result%expected_return = portfolio_mean(weights,mu)
    result%risk = portfolio_stddev(weights,sigma)
    if (result%risk > options%tolerance) then
      result%sharpe = (result%expected_return-options%risk_free)/result%risk
    end if
    call portfolio_returns(returns,weights,rp)
    result%expected_shortfall = historical_es(rp,options%alpha)
    if (allocated(c%initial_weights)) result%turnover = turnover(weights,c%initial_weights)
    result%hhi = hhi(weights)
    result%diversification = diversification(weights)
    call risk_contributions(weights,sigma,result%risk_contribution)
    result%feasible = check_constraints(weights,c,mu)
    result%converged = converged .and. result%feasible
    result%iterations = iters
    result%evaluations = evals
    if (result%feasible) then
      result%status = pa_success
      result%message = 'optimization completed'
      if (.not. converged) then
        result%status = pa_max_iterations
        result%message = 'feasible solution found before iteration limit'
      end if
    else
      result%status = pa_infeasible
      result%message = 'no feasible solution found'
    end if
  end subroutine fill_result

end module pa_optimizers
