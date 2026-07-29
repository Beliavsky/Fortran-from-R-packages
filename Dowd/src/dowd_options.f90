! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_options
  use dowd_kinds, only: dp
  use dowd_math, only: normal_cdf, normal_quantile, random_normal, sort_in_place
  use dowd_risk, only: historical_es
  implicit none
  private

  public :: black_scholes_call_price, black_scholes_put_price
  public :: long_call_var, short_call_var, long_put_var, short_put_var
  public :: black_scholes_call_es_sim, black_scholes_put_es_sim
  public :: american_put_price_binomial
  public :: american_put_var_binomial, american_put_es_binomial
  public :: american_put_es_sim

contains

  pure real(dp) function black_scholes_call_price(stock_price, strike, rate, sigma, maturity) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, sigma, maturity
    real(dp) :: d1, d2
    if (stock_price < 0.0_dp .or. strike < 0.0_dp .or. sigma < 0.0_dp .or. maturity < 0.0_dp) &
      error stop "black_scholes_call_price: invalid input"
    if (maturity <= tiny(1.0_dp) .or. sigma <= tiny(1.0_dp)) then
      value = max(stock_price-strike*exp(-rate*maturity),0.0_dp)
      return
    end if
    d1 = (log(stock_price/strike)+(rate+0.5_dp*sigma*sigma)*maturity)/(sigma*sqrt(maturity))
    d2 = d1-sigma*sqrt(maturity)
    value = stock_price*normal_cdf(d1)-strike*exp(-rate*maturity)*normal_cdf(d2)
  end function black_scholes_call_price

  pure real(dp) function black_scholes_put_price(stock_price, strike, rate, sigma, maturity) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, sigma, maturity
    real(dp) :: d1, d2
    if (stock_price < 0.0_dp .or. strike < 0.0_dp .or. sigma < 0.0_dp .or. maturity < 0.0_dp) &
      error stop "black_scholes_put_price: invalid input"
    if (maturity <= tiny(1.0_dp) .or. sigma <= tiny(1.0_dp)) then
      value = max(strike*exp(-rate*maturity)-stock_price,0.0_dp)
      return
    end if
    d1 = (log(stock_price/strike)+(rate+0.5_dp*sigma*sigma)*maturity)/(sigma*sqrt(maturity))
    d2 = d1-sigma*sqrt(maturity)
    value = strike*exp(-rate*maturity)*normal_cdf(-d2)-stock_price*normal_cdf(-d1)
  end function black_scholes_put_price

  pure real(dp) function long_call_var(stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days
    real(dp) :: t, hp, current_price, critical_stock, future_price
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (hp > t) error stop "long_call_var: holding period exceeds maturity"
    current_price = black_scholes_call_price(stock_price,strike,rate,sigma,t)
    critical_stock = exp(log(stock_price)+(drift-0.5_dp*sigma*sigma)*hp - &
                         normal_quantile(cl)*sigma*sqrt(hp))
    future_price = black_scholes_call_price(critical_stock,strike,rate,sigma,t-hp)
    value = current_price-future_price
  end function long_call_var

  pure real(dp) function short_call_var(stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days
    real(dp) :: t, hp, current_price, critical_stock, future_price
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (hp > t) error stop "short_call_var: holding period exceeds maturity"
    current_price = black_scholes_call_price(stock_price,strike,rate,sigma,t)
    critical_stock = exp(log(stock_price)+(drift-0.5_dp*sigma*sigma)*hp + &
                         normal_quantile(cl)*sigma*sqrt(hp))
    future_price = black_scholes_call_price(critical_stock,strike,rate,sigma,t-hp)
    value = future_price-current_price
  end function short_call_var

  pure real(dp) function long_put_var(stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days
    real(dp) :: t, hp, current_price, critical_stock, future_price
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (hp > t) error stop "long_put_var: holding period exceeds maturity"
    current_price = black_scholes_put_price(stock_price,strike,rate,sigma,t)
    critical_stock = exp(log(stock_price)+(drift-0.5_dp*sigma*sigma)*hp + &
                         normal_quantile(cl)*sigma*sqrt(hp))
    future_price = black_scholes_put_price(critical_stock,strike,rate,sigma,t-hp)
    value = current_price-future_price
  end function long_put_var

  pure real(dp) function short_put_var(stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, drift, sigma, maturity_days, cl, hp_days
    real(dp) :: t, hp, current_price, critical_stock, future_price
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (hp > t) error stop "short_put_var: holding period exceeds maturity"
    current_price = black_scholes_put_price(stock_price,strike,rate,sigma,t)
    critical_stock = exp(log(stock_price)+(drift-0.5_dp*sigma*sigma)*hp - &
                         normal_quantile(cl)*sigma*sqrt(hp))
    future_price = black_scholes_put_price(critical_stock,strike,rate,sigma,t-hp)
    value = future_price-current_price
  end function short_put_var

  real(dp) function black_scholes_call_es_sim(amount_invested, stock_price, strike, rate, drift, &
      sigma, maturity_days, number_trials, cl, hp_days) result(value)
    real(dp), intent(in) :: amount_invested, stock_price, strike, rate, drift, sigma
    real(dp), intent(in) :: maturity_days, cl, hp_days
    integer, intent(in) :: number_trials
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: t, hp, initial_price, number_options, new_stock, new_price, sign_position
    integer :: i
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (hp > t .or. number_trials <= 0) error stop "black_scholes_call_es_sim: invalid input"
    initial_price = black_scholes_call_price(stock_price,strike,rate,sigma,t)
    if (initial_price <= 0.0_dp) error stop "black_scholes_call_es_sim: zero initial option price"
    number_options = abs(amount_invested)/initial_price
    sign_position = merge(1.0_dp,-1.0_dp,amount_invested >= 0.0_dp)
    allocate(profit_loss(number_trials))
    do i = 1, number_trials
      new_stock = stock_price*exp((drift-0.5_dp*sigma*sigma)*hp+sigma*sqrt(hp)*random_normal())
      new_price = black_scholes_call_price(new_stock,strike,rate,sigma,t-hp)
      profit_loss(i) = sign_position*(new_price-initial_price)*number_options
    end do
    value = historical_es(profit_loss,cl)
  end function black_scholes_call_es_sim

  real(dp) function black_scholes_put_es_sim(amount_invested, stock_price, strike, rate, drift, &
      sigma, maturity_days, number_trials, cl, hp_days) result(value)
    real(dp), intent(in) :: amount_invested, stock_price, strike, rate, drift, sigma
    real(dp), intent(in) :: maturity_days, cl, hp_days
    integer, intent(in) :: number_trials
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: t, hp, initial_price, number_options, new_stock, new_price, sign_position
    integer :: i
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (hp > t .or. number_trials <= 0) error stop "black_scholes_put_es_sim: invalid input"
    initial_price = black_scholes_put_price(stock_price,strike,rate,sigma,t)
    if (initial_price <= 0.0_dp) error stop "black_scholes_put_es_sim: zero initial option price"
    number_options = abs(amount_invested)/initial_price
    sign_position = merge(1.0_dp,-1.0_dp,amount_invested >= 0.0_dp)
    allocate(profit_loss(number_trials))
    do i = 1, number_trials
      new_stock = stock_price*exp((drift-0.5_dp*sigma*sigma)*hp+sigma*sqrt(hp)*random_normal())
      new_price = black_scholes_put_price(new_stock,strike,rate,sigma,t-hp)
      profit_loss(i) = sign_position*(new_price-initial_price)*number_options
    end do
    value = historical_es(profit_loss,cl)
  end function black_scholes_put_es_sim

  real(dp) function american_put_price_binomial(stock_price, strike, rate, sigma, maturity_days, number_steps) result(value)
    real(dp), intent(in) :: stock_price, strike, rate, sigma, maturity_days
    integer, intent(in) :: number_steps
    real(dp), allocatable :: prices(:), option_values(:)
    real(dp) :: maturity, dt, u, d, growth, p, discount, exercise
    integer :: i, j, n
    n = number_steps
    if (n <= 0 .or. maturity_days < 0.0_dp) error stop "american_put_price_binomial: invalid input"
    maturity = maturity_days/360.0_dp
    if (maturity <= tiny(1.0_dp)) then
      value = max(strike-stock_price,0.0_dp)
      return
    end if
    dt = maturity/real(n,dp)
    u = exp(sigma*sqrt(dt))
    d = 1.0_dp/u
    growth = exp(rate*dt)
    p = (growth-d)/(u-d)
    p = min(max(p,0.0_dp),1.0_dp)
    discount = exp(-rate*dt)
    allocate(prices(n+1),option_values(n+1))
    do i = 0, n
      prices(i+1) = stock_price*u**i*d**(n-i)
      option_values(i+1) = max(strike-prices(i+1),0.0_dp)
    end do
    do j = n-1, 0, -1
      do i = 0, j
        prices(i+1) = stock_price*u**i*d**(j-i)
        option_values(i+1) = discount*(p*option_values(i+2)+(1.0_dp-p)*option_values(i+1))
        exercise = max(strike-prices(i+1),0.0_dp)
        option_values(i+1) = max(option_values(i+1),exercise)
      end do
    end do
    value = option_values(1)
  end function american_put_price_binomial

  subroutine american_put_holding_distribution(amount_invested, stock_price, strike, rate, sigma, &
      maturity_days, number_steps, hp_days, profit_loss, probabilities)
    real(dp), intent(in) :: amount_invested, stock_price, strike, rate, sigma, maturity_days, hp_days
    integer, intent(in) :: number_steps
    real(dp), allocatable, intent(out) :: profit_loss(:), probabilities(:)
    real(dp) :: maturity, hp, dt, u, d, growth, p, discount, initial_price
    real(dp) :: stock_node, number_options, sign_position
    integer :: n, m, i, remaining
    n = number_steps
    maturity = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    if (n <= 0 .or. hp < 0.0_dp .or. hp > maturity) error stop "american_put_holding_distribution: invalid input"
    dt = maturity/real(n,dp)
    m = min(n,max(0,nint(hp/dt)))
    u = exp(sigma*sqrt(dt))
    d = 1.0_dp/u
    growth = exp(rate*dt)
    p = min(max((growth-d)/(u-d),0.0_dp),1.0_dp)
    discount = exp(-rate*dt)
    initial_price = american_put_price_binomial(stock_price,strike,rate,sigma,maturity_days,n)
    if (initial_price <= 0.0_dp) error stop "american_put_holding_distribution: zero initial price"
    number_options = abs(amount_invested)/initial_price
    sign_position = merge(1.0_dp,-1.0_dp,amount_invested >= 0.0_dp)
    allocate(profit_loss(m+1),probabilities(m+1))
    do i = 0, m
      stock_node = stock_price*u**i*d**(m-i)
      remaining = n-m
      profit_loss(i+1) = sign_position*(american_put_price_binomial(stock_node,strike,rate,sigma, &
                              real(remaining,dp)*dt*360.0_dp,max(1,remaining))-initial_price)*number_options
      probabilities(i+1) = binomial_probability(m,i,p)
    end do
  contains
    pure real(dp) function binomial_probability(nn,kk,pp) result(prob)
      integer, intent(in) :: nn,kk
      real(dp), intent(in) :: pp
      integer :: l
      real(dp) :: coeff
      coeff = 1.0_dp
      do l = 1, min(kk,nn-kk)
        coeff = coeff*real(nn-l+1,dp)/real(l,dp)
      end do
      prob = coeff*pp**kk*(1.0_dp-pp)**(nn-kk)
    end function binomial_probability
  end subroutine american_put_holding_distribution

  real(dp) function weighted_loss_var(profit_loss, probabilities, cl) result(value)
    real(dp), intent(in) :: profit_loss(:), probabilities(:), cl
    real(dp), allocatable :: losses(:), probs(:)
    real(dp) :: cum, temp
    integer :: i, j, n
    n = size(profit_loss)
    allocate(losses(n),probs(n))
    losses = -profit_loss
    probs = probabilities/sum(probabilities)
    do i = 1, n-1
      do j = i+1, n
        if (losses(j) < losses(i)) then
          temp=losses(i); losses(i)=losses(j); losses(j)=temp
          temp=probs(i); probs(i)=probs(j); probs(j)=temp
        end if
      end do
    end do
    cum = 0.0_dp
    value = losses(n)
    do i = 1, n
      cum = cum+probs(i)
      if (cum >= cl) then
        value = losses(i)
        exit
      end if
    end do
  end function weighted_loss_var

  real(dp) function weighted_loss_es(profit_loss, probabilities, cl) result(value)
    real(dp), intent(in) :: profit_loss(:), probabilities(:), cl
    real(dp), allocatable :: losses(:), probs(:)
    real(dp) :: tail_prob, remaining, take, temp
    integer :: i, j, n
    n = size(profit_loss)
    allocate(losses(n),probs(n))
    losses = -profit_loss
    probs = probabilities/sum(probabilities)
    do i = 1, n-1
      do j = i+1, n
        if (losses(j) > losses(i)) then
          temp=losses(i); losses(i)=losses(j); losses(j)=temp
          temp=probs(i); probs(i)=probs(j); probs(j)=temp
        end if
      end do
    end do
    tail_prob = 1.0_dp-cl
    remaining = tail_prob
    value = 0.0_dp
    do i = 1, n
      take = min(remaining,probs(i))
      value = value+take*losses(i)
      remaining = remaining-take
      if (remaining <= 10.0_dp*epsilon(1.0_dp)) exit
    end do
    value = value/tail_prob
  end function weighted_loss_es

  real(dp) function american_put_var_binomial(amount_invested, stock_price, strike, rate, sigma, &
      maturity_days, number_steps, cl, hp_days) result(value)
    real(dp), intent(in) :: amount_invested, stock_price, strike, rate, sigma, maturity_days, cl, hp_days
    integer, intent(in) :: number_steps
    real(dp), allocatable :: profit_loss(:), probabilities(:)
    call american_put_holding_distribution(amount_invested,stock_price,strike,rate,sigma, &
         maturity_days,number_steps,hp_days,profit_loss,probabilities)
    value = weighted_loss_var(profit_loss,probabilities,cl)
  end function american_put_var_binomial

  real(dp) function american_put_es_binomial(amount_invested, stock_price, strike, rate, sigma, &
      maturity_days, number_steps, cl, hp_days) result(value)
    real(dp), intent(in) :: amount_invested, stock_price, strike, rate, sigma, maturity_days, cl, hp_days
    integer, intent(in) :: number_steps
    real(dp), allocatable :: profit_loss(:), probabilities(:)
    call american_put_holding_distribution(amount_invested,stock_price,strike,rate,sigma, &
         maturity_days,number_steps,hp_days,profit_loss,probabilities)
    value = weighted_loss_es(profit_loss,probabilities,cl)
  end function american_put_es_binomial

  real(dp) function american_put_es_sim(amount_invested, stock_price, strike, rate, drift, sigma, &
      maturity_days, number_trials, pricing_steps, cl, hp_days) result(value)
    real(dp), intent(in) :: amount_invested, stock_price, strike, rate, drift, sigma
    real(dp), intent(in) :: maturity_days, cl, hp_days
    integer, intent(in) :: number_trials, pricing_steps
    real(dp), allocatable :: profit_loss(:)
    real(dp) :: t, hp, initial_price, new_stock, new_price, number_options, sign_position
    integer :: i
    t = maturity_days/360.0_dp
    hp = hp_days/360.0_dp
    initial_price = american_put_price_binomial(stock_price,strike,rate,sigma,maturity_days,pricing_steps)
    number_options = abs(amount_invested)/initial_price
    sign_position = merge(1.0_dp,-1.0_dp,amount_invested >= 0.0_dp)
    allocate(profit_loss(number_trials))
    do i = 1, number_trials
      new_stock = stock_price*exp((drift-0.5_dp*sigma*sigma)*hp+sigma*sqrt(hp)*random_normal())
      new_price = american_put_price_binomial(new_stock,strike,rate,sigma,(t-hp)*360.0_dp,pricing_steps)
      profit_loss(i) = sign_position*(new_price-initial_price)*number_options
    end do
    value = historical_es(profit_loss,cl)
  end function american_put_es_sim

end module dowd_options
