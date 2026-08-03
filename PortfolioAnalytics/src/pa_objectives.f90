! SPDX-License-Identifier: GPL-3.0-only
module pa_objectives
  use pa_kinds, only: dp, pa_huge
  use pa_types, only: portfolio_constraints, portfolio_options, &
       obj_min_variance, obj_max_return, obj_quadratic_utility, obj_max_sharpe, &
       obj_min_es, obj_max_starr, obj_risk_parity, obj_min_semideviation, &
       obj_min_drawdown, obj_min_concentration, obj_min_csm, obj_max_csm_ratio, obj_min_eqs
  use pa_statistics, only: portfolio_returns, portfolio_mean, portfolio_variance, &
       portfolio_stddev, historical_es, semideviation, maximum_drawdown, turnover, &
       hhi, risk_contributions, conditional_second_moment, expected_quadratic_shortfall
  use pa_constraints, only: constraint_violation
  implicit none
  private
  public :: evaluate_portfolio_objective, objective_gradient_numeric

contains

  real(dp) function evaluate_portfolio_objective(weights, returns, mu, sigma, c, options, budgets) result(value)
    real(dp), intent(in) :: weights(:), returns(:,:), mu(:), sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: rp(:), rc(:), b(:)
    real(dp) :: mean_value, sd_value, risk_value, base, penalty, total_rc
    integer :: n

    n = size(weights)
    allocate(rp(size(returns,1)), rc(n), b(n))
    call portfolio_returns(returns,weights,rp)
    mean_value = portfolio_mean(weights,mu)
    sd_value = portfolio_stddev(weights,sigma)
    base = pa_huge
    select case(options%objective)
    case(obj_min_variance)
      base = portfolio_variance(weights,sigma)
    case(obj_max_return)
      base = -mean_value
    case(obj_quadratic_utility)
      base = 0.5_dp*options%risk_aversion*portfolio_variance(weights,sigma)-mean_value
    case(obj_max_sharpe)
      if (sd_value > options%tolerance) then
        base = -(mean_value-options%risk_free)/sd_value
      else
        base = pa_huge
      end if
    case(obj_min_es)
      base = historical_es(rp,options%alpha)
    case(obj_max_starr)
      risk_value = historical_es(rp,options%alpha)
      if (risk_value > options%tolerance) then
        base = -(mean_value-options%target_return)/risk_value
      else
        base = pa_huge
      end if
    case(obj_risk_parity)
      call risk_contributions(weights,sigma,rc)
      total_rc = sum(rc)
      if (present(budgets)) then
        if (size(budgets) == n .and. sum(budgets) > 0.0_dp) then
          b = budgets/sum(budgets)
        else
          b = 1.0_dp/real(n,dp)
        end if
      else
        b = 1.0_dp/real(n,dp)
      end if
      if (abs(total_rc) > options%tolerance) then
        base = sum((rc/total_rc-b)**2)
      else
        base = pa_huge
      end if
    case(obj_min_semideviation)
      base = semideviation(rp,options%target_return)
    case(obj_min_drawdown)
      base = maximum_drawdown(rp)
    case(obj_min_concentration)
      base = hhi(weights)
    case(obj_min_csm)
      base = conditional_second_moment(rp,options%alpha)
    case(obj_max_csm_ratio)
      risk_value = conditional_second_moment(rp,options%alpha)
      if (risk_value > options%tolerance) then
        base = -(mean_value-options%target_return)/risk_value
      else
        base = pa_huge
      end if
    case(obj_min_eqs)
      base = expected_quadratic_shortfall(rp,options%alpha)
    case default
      base = pa_huge
    end select

    if (allocated(c%initial_weights) .and. abs(options%turnover_aversion) > epsilon(1.0_dp)) then
      base = base + options%turnover_aversion*turnover(weights,c%initial_weights)
    end if
    if (abs(options%concentration_aversion) > epsilon(1.0_dp)) then
      base = base + options%concentration_aversion*hhi(weights)
    end if
    penalty = constraint_violation(weights,c,mu)
    value = base + options%penalty_scale*penalty
    if (.not. (value < pa_huge)) value = pa_huge
  end function evaluate_portfolio_objective

  subroutine objective_gradient_numeric(weights, returns, mu, sigma, c, options, gradient, budgets)
    real(dp), intent(in) :: weights(:), returns(:,:), mu(:), sigma(:,:)
    type(portfolio_constraints), intent(in) :: c
    type(portfolio_options), intent(in) :: options
    real(dp), intent(out) :: gradient(:)
    real(dp), intent(in), optional :: budgets(:)
    real(dp), allocatable :: wp(:), wm(:)
    real(dp) :: h, fp, fm
    integer :: i, n
    n = size(weights)
    allocate(wp(n),wm(n))
    do i = 1, n
      h = sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(weights(i)))
      wp = weights
      wm = weights
      wp(i) = wp(i)+h
      wm(i) = wm(i)-h
      if (present(budgets)) then
        fp = evaluate_portfolio_objective(wp,returns,mu,sigma,c,options,budgets)
        fm = evaluate_portfolio_objective(wm,returns,mu,sigma,c,options,budgets)
      else
        fp = evaluate_portfolio_objective(wp,returns,mu,sigma,c,options)
        fm = evaluate_portfolio_objective(wm,returns,mu,sigma,c,options)
      end if
      gradient(i) = (fp-fm)/(2.0_dp*h)
    end do
  end subroutine objective_gradient_numeric

end module pa_objectives
