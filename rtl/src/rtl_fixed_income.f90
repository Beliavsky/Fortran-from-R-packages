! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_fixed_income
  use rtl_kinds, only: dp
  use rtl_types, only: bond_result, npv_result, irs_result, commodity_weight_result
  use rtl_stats, only: natural_cubic_interpolate
  use rtl_calendar, only: add_months, swap_fut_weight
  implicit none
  private

  public :: bond_value, npv_value, interest_rate_swap, commodity_swap_prices
  public :: commodity_swap_from_calendar
  public :: bond, npv, swapIRS, swapCOM

  interface bond
    module procedure bond_value
  end interface bond

  interface npv
    module procedure npv_value
  end interface npv

  interface swapIRS
    module procedure interest_rate_swap
  end interface swapIRS

  interface swapCOM
    module procedure commodity_swap_prices
    module procedure commodity_swap_from_calendar
  end interface swapCOM

contains

  function bond_value(ytm, coupon_rate, maturity, payments_per_year) result(output)
    real(dp), intent(in) :: ytm, coupon_rate, maturity
    integer, intent(in) :: payments_per_year
    type(bond_result) :: output
    integer :: periods, i
    real(dp), parameter :: face_value = 100.0_dp

    if (payments_per_year < 1 .or. maturity <= 0.0_dp .or. &
        abs(maturity * real(payments_per_year, dp) - &
        real(nint(maturity * real(payments_per_year, dp)), dp)) > 1.0e-10_dp) then
      output%status%ok = .false.
      output%status%message = "maturity times payments_per_year must be a positive integer"
      return
    end if
    if (1.0_dp + ytm / real(payments_per_year, dp) <= 0.0_dp) then
      output%status%ok = .false.
      output%status%message = "yield produces invalid discount factors"
      return
    end if
    periods = nint(maturity * real(payments_per_year, dp))
    allocate(output%time(periods), output%cash_flow(periods), output%discount_factor(periods))
    allocate(output%present_value(periods), output%duration_contribution(periods))
    output%time = 0.0_dp
    output%cash_flow = 0.0_dp
    output%discount_factor = 0.0_dp
    output%present_value = 0.0_dp
    output%duration_contribution = 0.0_dp
    do i = 1, periods
      output%time(i) = real(i, dp) / real(payments_per_year, dp)
      if (i == periods) then
        output%cash_flow(i) = coupon_rate * face_value / real(payments_per_year, dp) + face_value
      else
        output%cash_flow(i) = coupon_rate * face_value / real(payments_per_year, dp)
      end if
      output%discount_factor(i) = 1.0_dp / &
        (1.0_dp + ytm / real(payments_per_year, dp))**i
    end do
    output%present_value = output%cash_flow * output%discount_factor
    output%price = sum(output%present_value)
    if (abs(output%price) > epsilon(1.0_dp)) then
      output%duration_contribution = output%present_value * output%time / output%price
      output%duration = sum(output%duration_contribution)
    end if
  end function bond_value

  function npv_value(initial_cost, periodic_cash_flow, cash_flow_frequency, terminal_value, &
      maturity, curve_times, curve_discounts, break_even, break_even_yield, &
      terminal_replaces_cash_flow) result(output)
    real(dp), intent(in) :: initial_cost, periodic_cash_flow, cash_flow_frequency
    real(dp), intent(in) :: terminal_value, maturity
    real(dp), intent(in) :: curve_times(:), curve_discounts(:)
    logical, intent(in), optional :: break_even, terminal_replaces_cash_flow
    real(dp), intent(in), optional :: break_even_yield
    type(npv_result) :: output
    integer :: periods, i
    logical :: flat_curve, replace_terminal
    real(dp) :: flat_yield

    if (cash_flow_frequency <= 0.0_dp .or. maturity <= 0.0_dp .or. &
        size(curve_times) /= size(curve_discounts) .or. size(curve_times) < 2) then
      output%status%ok = .false.
      output%status%message = "invalid NPV inputs"
      return
    end if
    periods = nint(maturity / cash_flow_frequency)
    if (abs(real(periods, dp) * cash_flow_frequency - maturity) > 1.0e-10_dp) then
      output%status%ok = .false.
      output%status%message = "maturity/cash_flow_frequency must be an integer"
      return
    end if
    flat_curve = .false.
    if (present(break_even)) flat_curve = break_even
    flat_yield = 0.0_dp
    if (present(break_even_yield)) flat_yield = break_even_yield
    replace_terminal = .true.
    if (present(terminal_replaces_cash_flow)) replace_terminal = terminal_replaces_cash_flow
    allocate(output%time(0:periods), output%cash_flow(0:periods))
    allocate(output%discount_factor(0:periods), output%present_value(0:periods))
    do i = 0, periods
      output%time(i) = real(i, dp) * cash_flow_frequency
      output%cash_flow(i) = periodic_cash_flow
      if (flat_curve) then
        output%discount_factor(i) = exp(-flat_yield * output%time(i))
      else
        output%discount_factor(i) = natural_cubic_interpolate(curve_times, curve_discounts, output%time(i))
      end if
    end do
    output%cash_flow(0) = initial_cost
    if (replace_terminal) then
      output%cash_flow(periods) = terminal_value
    else
      output%cash_flow(periods) = output%cash_flow(periods) + terminal_value
    end if
    output%present_value = output%cash_flow * output%discount_factor
    output%value = sum(output%present_value)
  end function npv_value

  function interest_rate_swap(trade_date, effective_date, maturity_date, notional, &
      pay_receive, fixed_rate, float_curve_times, float_curve_discounts, reset_months, &
      discount_curve_times, discount_curve_discounts, days_in_year) result(output)
    integer, intent(in) :: trade_date, effective_date, maturity_date, reset_months, days_in_year
    real(dp), intent(in) :: notional, fixed_rate
    character(len=*), intent(in) :: pay_receive
    real(dp), intent(in) :: float_curve_times(:), float_curve_discounts(:)
    real(dp), intent(in) :: discount_curve_times(:), discount_curve_discounts(:)
    type(irs_result) :: output
    integer, allocatable :: dates(:)
    real(dp), allocatable :: float_df(:), discount(:)
    integer :: n_payments, current_date, i, start_date, end_date, day_count
    real(dp) :: receiver_pv
    character(len=:), allocatable :: direction

    if (maturity_date <= effective_date .or. effective_date < trade_date .or. &
        notional < 0.0_dp .or. reset_months < 1 .or. &
        (days_in_year /= 360 .and. days_in_year /= 365)) then
      output%status%ok = .false.
      output%status%message = "invalid interest-rate-swap inputs"
      return
    end if
    if (size(float_curve_times) /= size(float_curve_discounts) .or. &
        size(discount_curve_times) /= size(discount_curve_discounts)) then
      output%status%ok = .false.
      output%status%message = "curve dimensions do not match"
      return
    end if
    current_date = effective_date
    n_payments = 0
    do while (current_date < maturity_date)
      current_date = min(add_months(current_date, reset_months), maturity_date)
      n_payments = n_payments + 1
    end do
    allocate(dates(n_payments), float_df(0:n_payments), discount(n_payments))
    current_date = effective_date
    do i = 1, n_payments
      current_date = min(add_months(current_date, reset_months), maturity_date)
      dates(i) = current_date
    end do
    allocate(output%payment_dates(n_payments), output%time(n_payments))
    allocate(output%fixed_leg(n_payments), output%floating_leg(n_payments), output%net(n_payments))
    output%payment_dates = dates
    float_df(0) = natural_cubic_interpolate(float_curve_times, float_curve_discounts, &
      real(effective_date - trade_date, dp) / 365.0_dp)
    do i = 1, n_payments
      output%time(i) = real(dates(i) - trade_date, dp) / 365.0_dp
      discount(i) = natural_cubic_interpolate(discount_curve_times, discount_curve_discounts, output%time(i))
      float_df(i) = natural_cubic_interpolate(float_curve_times, float_curve_discounts, output%time(i))
      if (i == 1) then
        start_date = effective_date
      else
        start_date = dates(i - 1)
      end if
      end_date = dates(i)
      day_count = end_date - start_date
      output%fixed_leg(i) = notional * fixed_rate * real(day_count, dp) / real(days_in_year, dp) * discount(i)
      output%floating_leg(i) = (float_df(i - 1) / float_df(i) - 1.0_dp) * notional * discount(i)
      output%net(i) = output%fixed_leg(i) - output%floating_leg(i)
    end do
    receiver_pv = sum(output%net)
    if (abs(receiver_pv) > 1.0e-14_dp) then
      output%duration = sum(real(dates - trade_date, dp) / real(days_in_year, dp) * output%net) / receiver_pv
    end if
    direction = lowercase(trim(pay_receive))
    if (direction == "rec" .or. direction == "receive") then
      output%present_value = receiver_pv
    else if (direction == "pay") then
      output%present_value = -receiver_pv
    else
      output%status%ok = .false.
      output%status%message = "pay_receive must be Pay, Rec, or Receive"
    end if
  end function interest_rate_swap

  pure function commodity_swap_prices(first_future, second_future, first_weight) result(prices)
    real(dp), intent(in) :: first_future(:), second_future(:), first_weight
    real(dp) :: prices(min(size(first_future), size(second_future)))
    integer :: n
    n = size(prices)
    prices = first_weight * first_future(1:n) + (1.0_dp - first_weight) * second_future(1:n)
  end function commodity_swap_prices


  function commodity_swap_from_calendar(first_future, second_future, month_start, &
      expiry_date, holidays) result(prices)
    real(dp), intent(in) :: first_future(:), second_future(:)
    integer, intent(in) :: month_start, expiry_date
    integer, intent(in), optional :: holidays(:)
    real(dp), allocatable :: prices(:)
    type(commodity_weight_result) :: weight_result

    if (present(holidays)) then
      weight_result = swap_fut_weight(month_start, expiry_date, holidays)
    else
      weight_result = swap_fut_weight(month_start, expiry_date)
    end if
    prices = commodity_swap_prices(first_future, second_future, weight_result%first_weight)
  end function commodity_swap_from_calendar

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

end module rtl_fixed_income
