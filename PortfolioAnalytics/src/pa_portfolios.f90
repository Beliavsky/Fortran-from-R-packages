! SPDX-License-Identifier: GPL-3.0-only
module pa_portfolios
  use pa_kinds, only: dp
  use pa_types, only: portfolio_constraints, portfolio_options, portfolio_result, &
       frontier_result, rebalancing_result, obj_min_variance, obj_max_return, &
       opt_differential_evolution
  use pa_constraints, only: random_feasible_portfolio, repair_weights
  use pa_optimizers, only: optimize_portfolio, equal_weight_portfolio
  use pa_statistics, only: sample_moments, portfolio_mean, portfolio_stddev
  implicit none
  private
  public :: random_portfolios, random_grid_portfolios, create_efficient_frontier
  public :: optimize_rebalancing

contains

  subroutine random_portfolios(c, nportfolios, weights, nfound, seed, mu)
    type(portfolio_constraints), intent(in) :: c
    integer, intent(in) :: nportfolios
    real(dp), intent(out) :: weights(:,:)
    integer, intent(out) :: nfound
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: mu(:)
    integer :: i, base_seed
    logical :: ok
    base_seed = 12345
    if (present(seed)) base_seed = seed
    nfound = 0
    do i = 1, min(nportfolios,size(weights,2))
      if (present(mu)) then
        call random_feasible_portfolio(c,weights(:,i),ok,base_seed+i,mu,10000)
      else
        call random_feasible_portfolio(c,weights(:,i),ok,base_seed+i,max_attempts=10000)
      end if
      if (ok) then
        nfound = nfound+1
      else
        weights(:,i) = 0.0_dp
      end if
    end do
  end subroutine random_portfolios

  subroutine random_grid_portfolios(c, nportfolios, spacing, weights, nfound, seed, mu)
    type(portfolio_constraints), intent(in) :: c
    integer, intent(in) :: nportfolios
    real(dp), intent(in) :: spacing
    real(dp), intent(out) :: weights(:,:)
    integer, intent(out) :: nfound
    integer, intent(in), optional :: seed
    real(dp), intent(in), optional :: mu(:)
    real(dp), allocatable :: raw(:), repaired(:)
    logical :: ok
    integer :: i, base_seed
    base_seed = 12345
    if (present(seed)) base_seed = seed
    allocate(raw(size(weights,1)),repaired(size(weights,1)))
    nfound = 0
    do i = 1, min(nportfolios,size(weights,2))
      if (present(mu)) then
        call random_feasible_portfolio(c,raw,ok,base_seed+i,mu,10000)
      else
        call random_feasible_portfolio(c,raw,ok,base_seed+i,max_attempts=10000)
      end if
      if (.not. ok) cycle
      raw = spacing*nint(raw/spacing)
      call repair_weights(raw,c,repaired,ok)
      if (.not. ok) cycle
      weights(:,i) = repaired
      nfound = nfound+1
    end do
    if (nfound < size(weights,2)) weights(:,nfound+1:) = 0.0_dp
  end subroutine random_grid_portfolios

  subroutine create_efficient_frontier(returns,c,options,npoints,frontier)
    real(dp), intent(in) :: returns(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    integer, intent(in) :: npoints
    type(frontier_result), intent(out) :: frontier
    type(portfolio_options) :: opt_min, opt_max, opt_i
    type(portfolio_constraints) :: ci
    type(portfolio_result) :: rmin, rmax, ri
    real(dp), allocatable :: mu(:), sigma(:,:)
    real(dp) :: lower, upper, target
    integer :: n, i

    n = size(returns,2)
    allocate(mu(n),sigma(n,n))
    call sample_moments(returns,mu,sigma)
    opt_min = options
    opt_min%objective = obj_min_variance
    call optimize_portfolio(returns,c,opt_min,rmin,mu,sigma)
    opt_max = options
    opt_max%objective = obj_max_return
    call optimize_portfolio(returns,c,opt_max,rmax,mu,sigma)
    lower = rmin%expected_return
    upper = max(lower,rmax%expected_return)
    allocate(frontier%weights(n,npoints),frontier%target_return(npoints), &
             frontier%expected_return(npoints),frontier%risk(npoints), &
             frontier%feasible(npoints))
    do i = 1, npoints
      if (npoints == 1) then
        target = lower
      else
        target = lower+(upper-lower)*real(i-1,dp)/real(npoints-1,dp)
      end if
      ci = c
      ci%return_target = target
      opt_i = options
      opt_i%objective = obj_min_variance
      if (target > lower+1.0e-10_dp) opt_i%optimizer = opt_differential_evolution
      call optimize_portfolio(returns,ci,opt_i,ri,mu,sigma)
      frontier%weights(:,i) = ri%weights
      frontier%target_return(i) = target
      frontier%expected_return(i) = ri%expected_return
      frontier%risk(i) = ri%risk
      frontier%feasible(i) = ri%feasible
    end do
  end subroutine create_efficient_frontier

  subroutine optimize_rebalancing(returns,rebalance,window,c,options,cost_rates,result)
    real(dp), intent(in) :: returns(:,:)
    logical, intent(in) :: rebalance(:)
    integer, intent(in) :: window
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(in), optional :: cost_rates(:)
    type(rebalancing_result), intent(out) :: result
    type(portfolio_constraints) :: ct
    type(portfolio_result) :: fit
    real(dp), allocatable :: current(:), previous(:), rates(:)
    real(dp) :: cost
    integer :: nobs, nassets, t, first
    logical :: ok

    nobs = size(returns,1)
    nassets = size(returns,2)
    allocate(result%weights(nassets,nobs),result%portfolio_return(nobs), &
             result%wealth(nobs),result%turnover(nobs),result%transaction_cost(nobs), &
             current(nassets),previous(nassets),rates(nassets))
    rates = 0.0_dp
    if (present(cost_rates)) rates = cost_rates
    call equal_weight_portfolio(c,current,ok)
    if (.not. ok) current = 1.0_dp/real(nassets,dp)
    previous = current
    result%wealth = 1.0_dp
    do t = 1, nobs
      previous = current
      if (t > 1 .and. rebalance(t)) then
        first = max(1,t-window)
        if (t-first >= 2) then
          ct = c
          if (allocated(ct%initial_weights)) deallocate(ct%initial_weights)
          allocate(ct%initial_weights(nassets))
          ct%initial_weights = previous
          call optimize_portfolio(returns(first:t-1,:),ct,options,fit)
          if (fit%feasible) current = fit%weights
        end if
      end if
      cost = sum(rates*abs(current-previous))
      result%transaction_cost(t) = cost
      result%turnover(t) = sum(abs(current-previous))/real(nassets,dp)
      result%weights(:,t) = current
      result%portfolio_return(t) = dot_product(current,returns(t,:))-cost
      if (t == 1) then
        result%wealth(t) = 1.0_dp+result%portfolio_return(t)
      else
        result%wealth(t) = result%wealth(t-1)*(1.0_dp+result%portfolio_return(t))
      end if
    end do
  end subroutine optimize_rebalancing

end module pa_portfolios
