! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_finance
  use r4gpf_kinds, only: dp
  use r4gpf_status, only: r4gpf_success, r4gpf_invalid_argument, r4gpf_dimension_error
  implicit none
  private
  public :: tax_assumptions, tax_result
  public :: purchasing_power, optimal_risky_asset_allocation, risk_adjusted_return
  public :: present_value_stream, utility, inverse_utility, certainty_equivalent_return
  public :: effective_tax_rates

  type :: tax_assumptions
    real(dp), allocatable :: turnover(:)
    real(dp), allocatable :: income_qualified(:)
    real(dp), allocatable :: capital_gains_long_term(:)
    real(dp), allocatable :: income(:)
    real(dp), allocatable :: capital_gains(:)
    real(dp), allocatable :: cost_basis(:)
  end type tax_assumptions

  type :: tax_result
    real(dp), allocatable :: blended_tax_rate_income(:)
    real(dp), allocatable :: blended_tax_rate_capital_gains(:)
    real(dp), allocatable :: preliquidation_aftertax_expected_return(:)
    real(dp), allocatable :: preliquidation_value(:)
    real(dp), allocatable :: capital_gain_taxed(:)
    real(dp), allocatable :: capital_gain_tax_paid(:)
    real(dp), allocatable :: postliquidation_value(:)
    real(dp), allocatable :: postliquidation_aftertax_expected_return(:)
    real(dp), allocatable :: effective_tax_rate(:)
    real(dp), allocatable :: aftertax_standard_deviation(:)
    integer :: status = r4gpf_success
  end type tax_result
contains

  elemental real(dp) function purchasing_power(x, years, real_interest_rate) result(value)
    real(dp), intent(in) :: x, years, real_interest_rate
    if (real_interest_rate < 0.0_dp) then
      value = x / (1.0_dp + abs(real_interest_rate))**years
    else
      value = x * (1.0_dp + real_interest_rate)**years
    end if
  end function purchasing_power

  elemental real(dp) function optimal_risky_asset_allocation(risky_return_mean, risky_return_sd, safe_return, &
      risk_aversion) result(allocation)
    real(dp), intent(in) :: risky_return_mean, risky_return_sd, safe_return, risk_aversion
    real(dp) :: denominator
    denominator = risk_aversion * risky_return_sd**2
    if (abs(denominator) <= tiny(1.0_dp)) then
      allocation = 0.0_dp
    else
      allocation = (risky_return_mean - safe_return) / denominator
      if (.not. (allocation > -huge(1.0_dp) .and. allocation < huge(1.0_dp))) allocation = 0.0_dp
    end if
  end function optimal_risky_asset_allocation

  elemental real(dp) function risk_adjusted_return(safe_return, risky_return_mean, risky_allocation, &
      risky_return_sd, risk_aversion) result(value)
    real(dp), intent(in) :: safe_return, risky_return_mean, risky_allocation
    real(dp), intent(in), optional :: risky_return_sd, risk_aversion
    real(dp) :: excess
    excess = risky_return_mean - safe_return
    if (.not. present(risky_return_sd) .and. .not. present(risk_aversion)) then
      value = safe_return + 0.5_dp * risky_allocation * excess
    else if (present(risky_return_sd) .and. present(risk_aversion)) then
      value = safe_return + risky_allocation * &
        (excess - 0.5_dp * risky_allocation * risk_aversion * risky_return_sd**2)
    else
      value = huge(1.0_dp)
    end if
  end function risk_adjusted_return

  subroutine present_value_stream(cashflow, discount_rate, present_values, status)
    real(dp), intent(in) :: cashflow(:), discount_rate
    real(dp), allocatable, intent(out) :: present_values(:)
    integer, intent(out), optional :: status
    integer :: i, j, n
    real(dp) :: discount

    n = size(cashflow)
    allocate(present_values(n))
    if (1.0_dp + discount_rate <= 0.0_dp) then
      present_values = huge(1.0_dp)
      if (present(status)) status = r4gpf_invalid_argument
      return
    end if
    do i = 1, n
      present_values(i) = 0.0_dp
      discount = 1.0_dp
      do j = i, n
        if (j > i) discount = discount * (1.0_dp + discount_rate)
        present_values(i) = present_values(i) + cashflow(j) / discount
      end do
    end do
    if (present(status)) status = r4gpf_success
  end subroutine present_value_stream

  elemental real(dp) function utility(x, parameter) result(value)
    real(dp), intent(in) :: x, parameter
    real(dp) :: xp
    xp = max(0.0_dp, x)
    if (xp <= tiny(1.0_dp) .or. abs(parameter) <= tiny(1.0_dp)) then
      value = 0.0_dp
    else if (abs(parameter - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) then
      value = log(xp)
    else
      value = (parameter / (parameter - 1.0_dp)) * &
        (xp**((parameter - 1.0_dp) / parameter) - 1.0_dp)
    end if
    if (.not. (value > -huge(1.0_dp) .and. value < huge(1.0_dp))) value = 0.0_dp
  end function utility

  elemental real(dp) function inverse_utility(value, parameter) result(x)
    real(dp), intent(in) :: value, parameter
    real(dp) :: base
    if (abs(parameter) <= tiny(1.0_dp)) then
      x = 1.0_dp
    else if (abs(parameter - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) then
      x = exp(value)
    else
      base = ((parameter - 1.0_dp) / parameter) * value + 1.0_dp
      if (base < 0.0_dp) then
        x = 0.0_dp
      else
        x = base**(parameter / (parameter - 1.0_dp))
      end if
    end if
  end function inverse_utility

  elemental real(dp) function certainty_equivalent_return(expected_return, variance, risk_tolerance) result(value)
    real(dp), intent(in) :: expected_return, variance, risk_tolerance
    if (risk_tolerance <= 0.0_dp) then
      value = -1.0_dp
    else
      value = (1.0_dp + expected_return) * exp(-variance / (2.0_dp * risk_tolerance)) - 1.0_dp
    end if
  end function certainty_equivalent_return

  subroutine effective_tax_rates(expected_return, standard_deviation, assumptions, tax_rate_ltcg, &
      tax_rate_ordinary_income, result, initial_value, investment_years)
    real(dp), intent(in) :: expected_return(:), standard_deviation(:)
    type(tax_assumptions), intent(in) :: assumptions
    real(dp), intent(in) :: tax_rate_ltcg, tax_rate_ordinary_income
    type(tax_result), intent(out) :: result
    real(dp), intent(in), optional :: initial_value, investment_years
    real(dp) :: iv, years
    integer :: n

    n = size(expected_return)
    if (size(standard_deviation) /= n .or. .not. assumptions_valid(assumptions, n)) then
      result%status = r4gpf_dimension_error
      return
    end if
    iv = 1000.0_dp
    years = 20.0_dp
    if (present(initial_value)) iv = initial_value
    if (present(investment_years)) years = investment_years
    allocate(result%blended_tax_rate_income(n), result%blended_tax_rate_capital_gains(n), &
      result%preliquidation_aftertax_expected_return(n), result%preliquidation_value(n), &
      result%capital_gain_taxed(n), result%capital_gain_tax_paid(n), result%postliquidation_value(n), &
      result%postliquidation_aftertax_expected_return(n), result%effective_tax_rate(n), &
      result%aftertax_standard_deviation(n))
    result%blended_tax_rate_income = assumptions%income_qualified * tax_rate_ltcg + &
      (1.0_dp - assumptions%income_qualified) * tax_rate_ordinary_income
    result%blended_tax_rate_capital_gains = assumptions%capital_gains_long_term * tax_rate_ltcg + &
      (1.0_dp - assumptions%capital_gains_long_term) * tax_rate_ordinary_income
    result%preliquidation_aftertax_expected_return = &
      (1.0_dp - result%blended_tax_rate_income) * assumptions%income + assumptions%capital_gains - &
      assumptions%turnover * (1.0_dp + assumptions%capital_gains - assumptions%cost_basis) * &
      result%blended_tax_rate_capital_gains
    result%preliquidation_value = iv * (1.0_dp + result%preliquidation_aftertax_expected_return)**years
    where (abs(expected_return) > tiny(1.0_dp))
      result%capital_gain_taxed = (assumptions%capital_gains / expected_return) * (1.0_dp - assumptions%turnover)
    elsewhere
      result%capital_gain_taxed = 0.0_dp
    end where
    result%capital_gain_tax_paid = (result%preliquidation_value - iv) * result%capital_gain_taxed * tax_rate_ltcg
    result%postliquidation_value = result%preliquidation_value - result%capital_gain_tax_paid
    result%postliquidation_aftertax_expected_return = &
      (result%postliquidation_value / iv)**(1.0_dp / years) - 1.0_dp
    where (abs(expected_return) > tiny(1.0_dp))
      result%effective_tax_rate = 1.0_dp - result%postliquidation_aftertax_expected_return / expected_return
    elsewhere
      result%effective_tax_rate = 0.0_dp
    end where
    result%aftertax_standard_deviation = (1.0_dp - result%effective_tax_rate) * standard_deviation
    result%status = r4gpf_success
  contains
    logical function assumptions_valid(a, expected_size) result(valid)
      type(tax_assumptions), intent(in) :: a
      integer, intent(in) :: expected_size
      valid = allocated(a%turnover) .and. allocated(a%income_qualified) .and. &
        allocated(a%capital_gains_long_term) .and. allocated(a%income) .and. &
        allocated(a%capital_gains) .and. allocated(a%cost_basis)
      if (.not. valid) return
      valid = size(a%turnover) == expected_size .and. size(a%income_qualified) == expected_size .and. &
        size(a%capital_gains_long_term) == expected_size .and. size(a%income) == expected_size .and. &
        size(a%capital_gains) == expected_size .and. size(a%cost_basis) == expected_size
    end function assumptions_valid
  end subroutine effective_tax_rates

end module r4gpf_finance
