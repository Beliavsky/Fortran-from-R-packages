! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_simulations
  use dowd_kinds, only: dp
  use dowd_math, only: mean_value, random_normal, random_lognormal, random_bernoulli
  use dowd_risk, only: historical_var, historical_es
  implicit none
  private

  public :: insurance_var_es
  public :: stop_loss_lognormal_var, filter_strategy_lognormal_var
  public :: default_risky_bond_var
  public :: dc_pension_var, db_pension_var
  public :: annuity_payment, present_value_annuity

contains

  subroutine insurance_var_es(log_mu, log_sigma, policies, claim_probability, premium_loading, &
      deductible, number_trials, cl, var_value, es_value)
    real(dp), intent(in) :: log_mu, log_sigma, claim_probability, premium_loading, deductible, cl
    integer, intent(in) :: policies, number_trials
    real(dp), intent(out) :: var_value, es_value
    real(dp), allocatable :: total_losses(:), profit_loss(:)
    real(dp) :: expected_loss
    integer :: i, j
    if (policies <= 0 .or. number_trials <= 0) error stop "insurance_var_es: invalid counts"
    allocate(total_losses(number_trials),profit_loss(number_trials))
    total_losses = 0.0_dp
    do j = 1, number_trials
      do i = 1, policies
        if (random_bernoulli(claim_probability) == 1) then
          total_losses(j) = total_losses(j)+max(random_lognormal(log_mu,log_sigma)-deductible,0.0_dp)
        end if
      end do
    end do
    expected_loss = mean_value(total_losses)
    profit_loss = -(total_losses-expected_loss-premium_loading*expected_loss/real(policies,dp))
    var_value = historical_var(profit_loss,cl)
    es_value = historical_es(profit_loss,cl)
  end subroutine insurance_var_es

  real(dp) function stop_loss_lognormal_var(mu, sigma, number_trials, loss_limit, cl, hp, number_steps) result(value)
    real(dp), intent(in) :: mu, sigma, loss_limit, cl, hp
    integer, intent(in) :: number_trials
    integer, intent(in), optional :: number_steps
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: dt, log_price, new_price, investment
    integer :: i, j, n
    n = 100
    if (present(number_steps)) n = max(1,number_steps)
    allocate(profit_loss(number_trials))
    do i = 1, number_trials
      log_price = 0.0_dp
      investment = 1.0_dp
      dt = hp/real(n,dp)
      do j = 1, n
        log_price = log_price+(mu-0.5_dp*sigma*sigma)*dt+sigma*sqrt(dt)*random_normal()
        new_price = exp(log_price)
        investment = max(new_price,1.0_dp-loss_limit)
      end do
      profit_loss(i) = investment-1.0_dp
    end do
    value = historical_var(profit_loss,cl)
  end function stop_loss_lognormal_var

  real(dp) function filter_strategy_lognormal_var(mu, sigma, number_trials, alpha, cl, hp, number_steps) result(value)
    real(dp), intent(in) :: mu, sigma, alpha, cl, hp
    integer, intent(in) :: number_trials
    integer, intent(in), optional :: number_steps
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: dt, log_price, new_price, equity_proportion, investment
    integer :: i, j, n
    n = 100
    if (present(number_steps)) n = max(1,number_steps)
    allocate(profit_loss(number_trials))
    dt = hp/real(n,dp)
    do i = 1, number_trials
      log_price = 0.0_dp
      investment = 1.0_dp
      do j = 1, n
        log_price = log_price+(mu-0.5_dp*sigma*sigma)*dt+sigma*sqrt(dt)*random_normal()
        new_price = exp(log_price)
        equity_proportion = 0.5_dp+alpha*(new_price-1.0_dp)
        investment = equity_proportion*new_price+(1.0_dp-equity_proportion)
      end do
      profit_loss(i) = investment-1.0_dp
    end do
    value = historical_var(profit_loss,cl)
  end function filter_strategy_lognormal_var

  real(dp) function default_risky_bond_var(spot_rate, risk_free_rate, coupon, rate_sigma, amount_invested, &
      recovery_rate, default_probability, number_trials, hp_days, cl) result(value)
    real(dp), intent(in) :: spot_rate, risk_free_rate, coupon, rate_sigma, amount_invested
    real(dp), intent(in) :: recovery_rate, default_probability, hp_days, cl
    integer, intent(in) :: number_trials
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: ann_hp, initial_value, number_bonds, random_rate, terminal_value, interim_value
    integer :: j, interim_default, terminal_default
    ann_hp = hp_days/360.0_dp
    initial_value = coupon*((1.0_dp+risk_free_rate)/(1.0_dp+spot_rate))**(ann_hp/2.0_dp) + &
                    (1.0_dp+coupon)/(1.0_dp+spot_rate)**ann_hp
    number_bonds = amount_invested/initial_value
    allocate(profit_loss(number_trials))
    do j = 1, number_trials
      random_rate = exp(log(max(spot_rate,tiny(1.0_dp)))+rate_sigma*sqrt(max(hp_days/2.0_dp,0.0_dp))*random_normal())
      interim_default = random_bernoulli(default_probability)
      terminal_default = random_bernoulli(default_probability)
      if (interim_default == 0) then
        if (terminal_default == 0) then
          terminal_value = coupon*(1.0_dp+risk_free_rate)**(ann_hp/2.0_dp)+(1.0_dp+coupon)
        else
          terminal_value = coupon*(1.0_dp+risk_free_rate)**(ann_hp/2.0_dp)+recovery_rate*(1.0_dp+coupon)
        end if
      else
        interim_value = recovery_rate*(coupon+(1.0_dp+coupon)/(1.0_dp+random_rate)**(ann_hp/2.0_dp))
        terminal_value = (1.0_dp+risk_free_rate)**(ann_hp/2.0_dp)*interim_value
      end if
      profit_loss(j) = number_bonds*(terminal_value-initial_value)
    end do
    value = historical_var(profit_loss,cl)
  end function default_risky_bond_var

  pure real(dp) function annuity_payment(rate, periods, present_value) result(payment)
    real(dp), intent(in) :: rate, present_value
    integer, intent(in) :: periods
    if (periods <= 0) error stop "annuity_payment: periods must be positive"
    if (abs(rate) < sqrt(epsilon(1.0_dp))) then
      payment = present_value/real(periods,dp)
    else
      payment = present_value*rate/(1.0_dp-(1.0_dp+rate)**(-periods))
    end if
  end function annuity_payment

  pure real(dp) function present_value_annuity(rate, periods, cashflow) result(value)
    real(dp), intent(in) :: rate, cashflow
    integer, intent(in) :: periods
    if (periods <= 0) error stop "present_value_annuity: periods must be positive"
    if (abs(rate) < sqrt(epsilon(1.0_dp))) then
      value = cashflow*real(periods,dp)
    else
      value = cashflow/rate*(1.0_dp-(1.0_dp+rate)**(-periods))
    end if
  end function present_value_annuity

  real(dp) function dc_pension_var(mu, sigma, unemployment_probability, life_expectancy, number_trials, cl) result(value)
    real(dp), intent(in) :: mu, sigma, unemployment_probability, life_expectancy, cl
    integer, intent(in) :: number_trials
    real(dp), allocatable :: pension_ratio(:)
    real(dp) :: fund, income, contribution, terminal_income, pension
    integer :: i, j, periods
    periods = max(1,int(life_expectancy)-65)
    allocate(pension_ratio(number_trials))
    terminal_income = (1.0_dp-unemployment_probability)*25.0_dp*exp(0.02_dp*39.0_dp)
    do j = 1, number_trials
      fund = 0.15_dp*25.0_dp
      do i = 2, 40
        income = 25.0_dp*exp(0.02_dp*real(i-1,dp))*real(random_bernoulli(1.0_dp-unemployment_probability),dp)
        contribution = 0.15_dp*income
        fund = contribution+fund*(1.0_dp+random_normal(mu,sigma))
      end do
      pension = annuity_payment(0.04_dp,periods,fund)
      pension_ratio(j) = pension/terminal_income
    end do
    value = -historical_var(pension_ratio,cl)
  end function dc_pension_var

  real(dp) function db_pension_var(mu, sigma, unemployment_probability, life_expectancy, number_trials, cl) result(value)
    real(dp), intent(in) :: mu, sigma, unemployment_probability, life_expectancy, cl
    integer, intent(in) :: number_trials
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: fund, income, contribution, terminal_income, pension, implied_fund
    integer :: i, j, years_contributed, periods
    periods = max(1,int(life_expectancy)-65)
    terminal_income = (1.0_dp-unemployment_probability)*25.0_dp*exp(0.02_dp*39.0_dp)
    allocate(profit_loss(number_trials))
    do j = 1, number_trials
      fund = 0.15_dp*25.0_dp
      years_contributed = 1
      do i = 2, 40
        if (random_bernoulli(1.0_dp-unemployment_probability) == 1) then
          income = 25.0_dp*exp(0.02_dp*real(i-1,dp))
          years_contributed = years_contributed+1
        else
          income = 0.0_dp
        end if
        contribution = 0.15_dp*income
        fund = contribution+fund*(1.0_dp+random_normal(mu,sigma))
      end do
      pension = real(years_contributed,dp)/40.0_dp*terminal_income
      implied_fund = present_value_annuity(0.04_dp,periods,pension)
      profit_loss(j) = fund-implied_fund
    end do
    value = historical_var(profit_loss,cl)
  end function db_pension_var

end module dowd_simulations
