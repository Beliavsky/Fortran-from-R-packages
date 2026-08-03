! SPDX-License-Identifier: MIT
module greeks_payoffs
  use greeks_kinds, only: dp
  implicit none
  private
  public :: payoff_value, payoff_derivative, valid_standard_payoff

contains

  pure logical function valid_standard_payoff(payoff) result(valid)
    character(len=*), intent(in) :: payoff
    valid = trim(payoff) == 'call' .or. trim(payoff) == 'put' .or. &
      trim(payoff) == 'digital_call' .or. trim(payoff) == 'digital_put' .or. &
      trim(payoff) == 'cash_or_nothing_call' .or. &
      trim(payoff) == 'cash_or_nothing_put' .or. &
      trim(payoff) == 'asset_or_nothing_call' .or. &
      trim(payoff) == 'asset_or_nothing_put'
  end function valid_standard_payoff

  pure real(dp) function payoff_value(x, strike, payoff) result(value)
    real(dp), intent(in) :: x, strike
    character(len=*), intent(in) :: payoff
    select case (trim(payoff))
    case ('call')
      value = max(0.0_dp, x - strike)
    case ('put')
      value = max(0.0_dp, strike - x)
    case ('digital_call', 'cash_or_nothing_call')
      value = merge(1.0_dp, 0.0_dp, x >= strike)
    case ('digital_put', 'cash_or_nothing_put')
      value = merge(1.0_dp, 0.0_dp, x <= strike)
    case ('asset_or_nothing_call')
      value = merge(x, 0.0_dp, x >= strike)
    case ('asset_or_nothing_put')
      value = merge(x, 0.0_dp, x <= strike)
    case default
      value = 0.0_dp
    end select
  end function payoff_value

  pure real(dp) function payoff_derivative(x, strike, payoff) result(value)
    real(dp), intent(in) :: x, strike
    character(len=*), intent(in) :: payoff
    select case (trim(payoff))
    case ('call')
      value = merge(1.0_dp, 0.0_dp, x > strike)
    case ('put')
      value = merge(-1.0_dp, 0.0_dp, x < strike)
    case default
      value = 0.0_dp
    end select
  end function payoff_derivative

end module greeks_payoffs
