! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_options
  use rtl_kinds, only: dp
  use rtl_types, only: option_result, spread_option_result, crr_tree_result
  use rtl_stats, only: normal_cdf, normal_pdf
  implicit none
  private

  public :: gbs_option, crr_option, crr_euro_tree
  public :: spread_option, barrier_spread_option
  public :: CRROption, CRReuro, GBSOption, spreadOption, barrierSpreadOption

  interface CRROption
    module procedure crr_option
  end interface CRROption

  interface CRReuro
    module procedure crr_euro_tree
  end interface CRReuro

  interface GBSOption
    module procedure gbs_option
  end interface GBSOption

  interface spreadOption
    module procedure spread_option
  end interface spreadOption

  interface barrierSpreadOption
    module procedure barrier_spread_option
  end interface barrierSpreadOption

contains

  function gbs_option(spot, strike, maturity, rate, carry, volatility, option_type) result(output)
    real(dp), intent(in) :: spot, strike, maturity, rate, carry, volatility
    character(len=*), intent(in) :: option_type
    type(option_result) :: output
    real(dp) :: sigma, d1, d2, discount_asset, discount_strike, sqrt_t
    character(len=:), allocatable :: kind_name

    kind_name = lowercase(trim(option_type))
    if (spot <= 0.0_dp .or. strike <= 0.0_dp .or. maturity < 0.0_dp .or. volatility < 0.0_dp) then
      call fail(output, "invalid generalized Black-Scholes inputs")
      return
    end if
    if (kind_name /= "call" .and. kind_name /= "put") then
      call fail(output, "option type must be call or put")
      return
    end if
    if (maturity <= epsilon(1.0_dp)) then
      if (kind_name == "call") then
        output%price = max(spot - strike, 0.0_dp)
        if (spot > strike) output%delta = 1.0_dp
      else
        output%price = max(strike - spot, 0.0_dp)
        if (spot < strike) output%delta = -1.0_dp
      end if
      return
    end if
    sigma = max(volatility, 1.0e-16_dp)
    sqrt_t = sqrt(maturity)
    d1 = (log(spot / strike) + (carry + 0.5_dp * sigma * sigma) * maturity) / (sigma * sqrt_t)
    d2 = d1 - sigma * sqrt_t
    discount_asset = exp((carry - rate) * maturity)
    discount_strike = exp(-rate * maturity)

    output%gamma = discount_asset * normal_pdf(d1) / (spot * sigma * sqrt_t)
    output%vega = spot * sqrt_t * discount_asset * normal_pdf(d1)
    if (kind_name == "call") then
      output%price = spot * discount_asset * normal_cdf(d1) - &
        strike * discount_strike * normal_cdf(d2)
      output%delta = discount_asset * normal_cdf(d1)
      output%theta = -spot * normal_pdf(d1) * sigma * discount_asset / (2.0_dp * sqrt_t) - &
        rate * strike * discount_strike * normal_cdf(d2) + &
        (carry - rate) * spot * discount_asset * normal_cdf(d1)
      output%rho = strike * maturity * discount_strike * normal_cdf(d2)
    else
      output%price = strike * discount_strike * normal_cdf(-d2) - &
        spot * discount_asset * normal_cdf(-d1)
      output%delta = -discount_asset * normal_cdf(-d1)
      output%theta = -spot * normal_pdf(d1) * sigma * discount_asset / (2.0_dp * sqrt_t) + &
        rate * strike * discount_strike * normal_cdf(-d2) - &
        (carry - rate) * spot * discount_asset * normal_cdf(-d1)
      output%rho = -strike * maturity * discount_strike * normal_cdf(-d2)
    end if
  end function gbs_option

  function crr_option(spot, strike, volatility, rate, carry, maturity, steps, &
      option_type, option_style) result(output)
    real(dp), intent(in) :: spot, strike, volatility, rate, carry, maturity
    integer, intent(in) :: steps
    character(len=*), intent(in) :: option_type, option_style
    type(option_result) :: output
    real(dp), allocatable :: values(:)
    real(dp) :: dt, up, down, probability, discount, hold_value, exercise_value
    integer :: i, j
    character(len=:), allocatable :: kind_name, style_name

    kind_name = lowercase(trim(option_type))
    style_name = lowercase(trim(option_style))
    if (steps < 1 .or. spot <= 0.0_dp .or. strike < 0.0_dp .or. maturity < 0.0_dp) then
      call fail(output, "invalid CRR inputs")
      return
    end if
    if (kind_name /= "call" .and. kind_name /= "put") then
      call fail(output, "option type must be call or put")
      return
    end if
    if (style_name /= "european" .and. style_name /= "american") then
      call fail(output, "option style must be european or american")
      return
    end if
    if (maturity <= epsilon(1.0_dp)) then
      output%price = payoff(spot, strike, kind_name)
      return
    end if

    dt = maturity / real(steps, dp)
    up = exp(max(volatility, 1.0e-16_dp) * sqrt(dt))
    down = 1.0_dp / up
    probability = (exp(carry * dt) - down) / (up - down)
    if (probability < 0.0_dp .or. probability > 1.0_dp) then
      output%status%message = "CRR risk-neutral probability lies outside [0,1]"
    end if
    discount = exp(-rate * dt)
    allocate(values(0:steps))
    do i = 0, steps
      values(i) = payoff(spot * up**i * down**(steps - i), strike, kind_name)
    end do
    do i = steps - 1, 0, -1
      do j = 0, i
        hold_value = discount * (probability * values(j + 1) + (1.0_dp - probability) * values(j))
        if (style_name == "american") then
          exercise_value = payoff(spot * up**j * down**(i - j), strike, kind_name)
          values(j) = max(hold_value, exercise_value)
        else
          values(j) = hold_value
        end if
      end do
    end do
    output%price = values(0)
  end function crr_option

  function crr_euro_tree(spot, strike, volatility, rate, maturity, steps, option_type) result(output)
    real(dp), intent(in) :: spot, strike, volatility, rate, maturity
    integer, intent(in) :: steps
    character(len=*), intent(in) :: option_type
    type(crr_tree_result) :: output
    real(dp) :: dt, up, down, probability, discount
    integer :: i, j
    character(len=:), allocatable :: kind_name

    kind_name = lowercase(trim(option_type))
    if (steps < 1 .or. maturity <= 0.0_dp .or. volatility <= 0.0_dp) then
      output%status%ok = .false.
      output%status%message = "invalid CRR tree inputs"
      return
    end if
    if (kind_name /= "call" .and. kind_name /= "put") then
      output%status%ok = .false.
      output%status%message = "option type must be call or put"
      return
    end if
    allocate(output%asset(0:steps, 0:steps), output%option(0:steps, 0:steps))
    output%asset = 0.0_dp
    output%option = 0.0_dp
    dt = maturity / real(steps, dp)
    up = exp(volatility * sqrt(dt))
    down = exp(-volatility * sqrt(dt))
    probability = (exp(rate * dt) - down) / (up - down)
    discount = exp(-rate * dt)
    do i = 0, steps
      do j = 0, i
        output%asset(i, j) = spot * up**j * down**(i - j)
      end do
    end do
    do j = 0, steps
      output%option(steps, j) = payoff(output%asset(steps, j), strike, kind_name)
    end do
    do i = steps - 1, 0, -1
      do j = 0, i
        output%option(i, j) = discount * ((1.0_dp - probability) * output%option(i + 1, j) + &
          probability * output%option(i + 1, j + 1))
      end do
    end do
    output%price = output%option(0, 0)
    if (volatility <= sqrt(dt) * rate) then
      output%status%message = "sigma < r*sqrt(dt); upstream warns against using this tree"
    end if
  end function crr_euro_tree

  function spread_option(f1, f2, strike, sigma1, sigma2, correlation, maturity, rate, &
      option_type) result(output)
    real(dp), intent(in) :: f1, f2, strike, sigma1, sigma2, correlation, maturity, rate
    character(len=*), intent(in) :: option_type
    type(spread_option_result) :: output
    real(dp), parameter :: epsilon = 1.0e-10_dp
    real(dp) :: w1, w2, sigma_eff, d1, d2, discount, sqrt_t
    character(len=:), allocatable :: kind_name

    kind_name = lowercase(trim(option_type))
    if (kind_name /= "call" .and. kind_name /= "put") then
      call fail_spread(output, "option type must be call or put")
      return
    end if
    if (maturity <= 0.0_dp .or. sigma1 < 0.0_dp .or. sigma2 < 0.0_dp .or. &
        abs(correlation) > 1.0_dp) then
      call fail_spread(output, "invalid spread-option inputs")
      return
    end if
    if (f2 + epsilon <= 0.0_dp .or. f1 + strike + epsilon <= 0.0_dp) then
      call fail_spread(output, "Kirk approximation requires positive ratio terms")
      return
    end if
    w1 = f1 / (f1 + f2 + epsilon)
    w2 = f2 / (f1 + f2 + epsilon)
    sigma_eff = sqrt(max(epsilon, sigma1**2 * w1**2 + sigma2**2 * w2**2 - &
      2.0_dp * correlation * sigma1 * sigma2 * w1 * w2 + epsilon))
    sqrt_t = sqrt(maturity)
    d1 = (log((f2 + epsilon) / (f1 + strike + epsilon)) + &
      0.5_dp * sigma_eff**2 * maturity) / (sigma_eff * sqrt(maturity + epsilon))
    d2 = d1 - sigma_eff * sqrt(maturity + epsilon)
    discount = exp(-rate * maturity)

    if (kind_name == "call") then
      output%price = discount * (f2 * normal_cdf(d1) - (f1 + strike) * normal_cdf(d2))
      output%delta_f1 = -discount * normal_cdf(d2)
      output%delta_f2 = discount * normal_cdf(d1)
    else
      output%price = discount * ((f1 + strike) * normal_cdf(-d2) - f2 * normal_cdf(-d1))
      output%delta_f1 = discount * normal_cdf(-d2)
      output%delta_f2 = -discount * normal_cdf(-d1)
    end if
    output%gamma_f1 = discount * normal_pdf(d2) / ((f1 + strike) * sigma_eff * sqrt_t)
    output%gamma_f2 = discount * normal_pdf(d1) / (f2 * sigma_eff * sqrt_t)
    output%gamma_cross = -output%gamma_f1
    output%vega_1 = discount * f2 * sqrt_t * normal_pdf(d1) * &
      (sigma1 * w1**2 - correlation * sigma2 * w1 * w2) / sigma_eff
    output%vega_2 = discount * f2 * sqrt_t * normal_pdf(d1) * &
      (sigma2 * w2**2 - correlation * sigma1 * w1 * w2) / sigma_eff
    if (kind_name == "call") then
      output%theta = -discount * (f2 * normal_pdf(d1) * sigma_eff / (2.0_dp * sqrt_t) - &
        (f1 + strike) * normal_pdf(d2) * sigma_eff / (2.0_dp * sqrt_t) + &
        rate * (f2 * normal_cdf(d1) - (f1 + strike) * normal_cdf(d2)))
    else
      output%theta = -discount * (f2 * normal_pdf(d1) * sigma_eff / (2.0_dp * sqrt_t) - &
        (f1 + strike) * normal_pdf(d2) * sigma_eff / (2.0_dp * sqrt_t) - &
        rate * ((f1 + strike) * normal_cdf(-d2) - f2 * normal_cdf(-d1)))
    end if
    output%rho = maturity * output%price
  end function spread_option

  function barrier_spread_option(f1, f2, strike, barrier, sigma1, sigma2, correlation, &
      maturity, rate, option_type, barrier_type, monitoring, upstream_round) result(output)
    real(dp), intent(in) :: f1, f2, strike, barrier, sigma1, sigma2, correlation, maturity, rate
    character(len=*), intent(in) :: option_type, barrier_type
    character(len=*), intent(in), optional :: monitoring
    logical, intent(in), optional :: upstream_round
    type(spread_option_result) :: output
    real(dp), parameter :: epsilon = 1.0e-10_dp
    real(dp) :: spread, sigma_spread, d1, d2, b1, b2, discount, power_value
    real(dp) :: vanilla, reflection, gamma_base, raw_price
    character(len=:), allocatable :: kind_name, barrier_name, monitor_name
    logical :: round_price

    kind_name = lowercase(trim(option_type))
    barrier_name = lowercase(trim(barrier_type))
    monitor_name = "continuous"
    if (present(monitoring)) monitor_name = lowercase(trim(monitoring))
    round_price = .false.
    if (present(upstream_round)) round_price = upstream_round

    if (sigma1 <= 0.0_dp .or. sigma2 <= 0.0_dp .or. abs(correlation) > 1.0_dp .or. &
        maturity < 0.0_dp) then
      call fail_spread(output, "invalid barrier-spread inputs")
      return
    end if
    if ((kind_name == "call" .and. barrier_name /= "uo") .or. &
        (kind_name == "put" .and. barrier_name /= "do")) then
      call fail_spread(output, "only up-and-out calls and down-and-out puts are supported")
      return
    end if
    if (monitor_name /= "continuous" .and. monitor_name /= "terminal") then
      call fail_spread(output, "monitoring must be continuous or terminal")
      return
    end if
    ! The upstream parameter is validated but is not used by its formula.
    output%monitoring_used = .false.
    spread = f2 - f1
    sigma_spread = sqrt(max(0.0_dp, sigma1**2 + sigma2**2 - &
      2.0_dp * correlation * sigma1 * sigma2))

    if ((kind_name == "call" .and. spread >= barrier) .or. &
        (kind_name == "put" .and. spread <= barrier)) return
    if (maturity < epsilon) then
      if (kind_name == "call") output%price = max(0.0_dp, spread - strike)
      if (kind_name == "put") output%price = max(0.0_dp, strike - spread)
      return
    end if
    if (spread <= 0.0_dp .or. strike <= 0.0_dp .or. barrier <= 0.0_dp .or. sigma_spread <= 0.0_dp) then
      call fail_spread(output, "barrier formula requires positive spread, strike, barrier, and spread volatility")
      return
    end if

    d1 = (log(spread / strike) + 0.5_dp * sigma_spread**2 * maturity) / &
      (sigma_spread * sqrt(maturity))
    d2 = d1 - sigma_spread * sqrt(maturity)
    b1 = (log(spread / barrier) + 0.5_dp * sigma_spread**2 * maturity) / &
      (sigma_spread * sqrt(maturity))
    b2 = b1 - sigma_spread * sqrt(maturity)
    discount = exp(-rate * maturity)

    if (kind_name == "call") then
      power_value = (barrier / spread)**2
      vanilla = discount * (spread * normal_cdf(d1) - strike * normal_cdf(d2))
      reflection = discount * power_value * (spread * normal_cdf(b1) - strike * normal_cdf(b2))
      raw_price = max(0.0_dp, vanilla - reflection)
      output%delta_f1 = -discount * (normal_cdf(d1) - power_value * normal_cdf(b1))
      output%delta_f2 = -output%delta_f1
      gamma_base = discount * (normal_pdf(d1) - power_value * normal_pdf(b1)) / &
        (spread * sigma_spread * sqrt(maturity))
      output%vega_1 = spread * discount * sqrt(maturity) * &
        (normal_pdf(d1) - power_value * normal_pdf(b1)) * &
        (sigma1 - correlation * sigma2) / sigma_spread
      output%vega_2 = spread * discount * sqrt(maturity) * &
        (normal_pdf(d1) - power_value * normal_pdf(b1)) * &
        (sigma2 - correlation * sigma1) / sigma_spread
      output%theta = -0.5_dp * sigma_spread**2 * spread**2 * gamma_base - &
        rate * (spread * normal_cdf(d1) - strike * normal_cdf(d2)) * discount + &
        rate * power_value * (spread * normal_cdf(b1) - strike * normal_cdf(b2)) * discount
      output%rho = maturity * strike * discount * &
        (normal_cdf(d2) - power_value * normal_cdf(b2))
    else
      power_value = (spread / barrier)**2
      vanilla = discount * (strike * normal_cdf(-d2) - spread * normal_cdf(-d1))
      reflection = discount * power_value * (strike * normal_cdf(-b2) - spread * normal_cdf(-b1))
      raw_price = max(0.0_dp, vanilla - reflection)
      output%delta_f1 = discount * (normal_cdf(-d1) - power_value * normal_cdf(-b1))
      output%delta_f2 = -output%delta_f1
      gamma_base = discount * (normal_pdf(d1) - power_value * normal_pdf(b1)) / &
        (spread * sigma_spread * sqrt(maturity))
      output%vega_1 = spread * discount * sqrt(maturity) * &
        (normal_pdf(d1) - power_value * normal_pdf(b1)) * &
        (sigma1 - correlation * sigma2) / sigma_spread
      output%vega_2 = spread * discount * sqrt(maturity) * &
        (normal_pdf(d1) - power_value * normal_pdf(b1)) * &
        (sigma2 - correlation * sigma1) / sigma_spread
      output%theta = -0.5_dp * sigma_spread**2 * spread**2 * gamma_base + &
        rate * (strike * normal_cdf(-d2) - spread * normal_cdf(-d1)) * discount - &
        rate * power_value * (strike * normal_cdf(-b2) - spread * normal_cdf(-b1)) * discount
      output%rho = -maturity * strike * discount * &
        (normal_cdf(-d2) - power_value * normal_cdf(-b2))
    end if
    output%gamma_f1 = gamma_base
    output%gamma_f2 = gamma_base
    output%gamma_cross = -gamma_base
    if (round_price) then
      output%price = real(nint(100.0_dp * raw_price), dp) / 100.0_dp
    else
      output%price = raw_price
    end if
  end function barrier_spread_option

  pure real(dp) function payoff(underlying, strike, kind_name) result(value)
    real(dp), intent(in) :: underlying, strike
    character(len=*), intent(in) :: kind_name
    if (kind_name == "call") then
      value = max(underlying - strike, 0.0_dp)
    else
      value = max(strike - underlying, 0.0_dp)
    end if
  end function payoff

  pure function lowercase(text) result(output)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: output
    integer :: i, code
    output = text
    do i = 1, len(text)
      code = iachar(output(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
    end do
  end function lowercase

  subroutine fail(output, message)
    type(option_result), intent(inout) :: output
    character(len=*), intent(in) :: message
    output%status%ok = .false.
    output%status%message = message
  end subroutine fail

  subroutine fail_spread(output, message)
    type(spread_option_result), intent(inout) :: output
    character(len=*), intent(in) :: message
    output%status%ok = .false.
    output%status%message = message
  end subroutine fail_spread

end module rtl_options
