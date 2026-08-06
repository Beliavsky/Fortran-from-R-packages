module trading_trades
  use trading_kinds, only : dp, str_len
  use trading_stats, only : normal_cdf
  implicit none
  private

  type, public :: trade_t
    character(len=str_len) :: class_name = "Trade"
    character(len=str_len) :: trade_group = ""
    character(len=str_len) :: trade_type = ""
    character(len=str_len) :: subclass = ""
    real(dp) :: notional = 0.0_dp
    real(dp) :: mtm = 0.0_dp
    character(len=str_len) :: currency = ""
    real(dp) :: si = 0.0_dp
    real(dp) :: ei = 0.0_dp
    character(len=str_len) :: buy_sell = ""
    character(len=str_len) :: external_id = ""
    character(len=str_len) :: counterparty = ""
    character(len=str_len) :: nickname = ""
    character(len=str_len) :: exotic_type = ""
    character(len=str_len) :: netting_set = ""
    character(len=str_len) :: underlying_instrument = ""
    logical :: simplified = .false.
    logical :: use_simplified = .false.
    character(len=str_len) :: ccy_pair = ""
    character(len=str_len) :: ccy_paying = ""
    real(dp) :: amount_paying = 0.0_dp
    character(len=str_len) :: ccy_receiving = ""
    real(dp) :: amount_receiving = 0.0_dp
    character(len=str_len) :: base_ccy = ""
    character(len=str_len) :: option_type = ""
    real(dp) :: underlying_price = 0.0_dp
    real(dp) :: strike_price = 0.0_dp
    character(len=str_len) :: delivery_type = ""
    real(dp) :: maturity = 0.0_dp
    character(len=str_len) :: exotic_option_type = ""
    character(len=str_len) :: pay_leg_type = ""
    character(len=str_len) :: pay_leg_ref = ""
    character(len=str_len) :: pay_leg_tenor = ""
    real(dp) :: pay_leg_rate = 0.0_dp
    character(len=str_len) :: rec_leg_type = ""
    character(len=str_len) :: rec_leg_ref = ""
    character(len=str_len) :: rec_leg_tenor = ""
    real(dp) :: rec_leg_rate = 0.0_dp
    real(dp) :: vol_strike = 0.0_dp
    real(dp) :: annualization_factor = 0.0_dp
    real(dp) :: vega_notional = 0.0_dp
    character(len=str_len) :: commodity_type = ""
    character(len=str_len) :: reference_entity = ""
    real(dp) :: cdo_attach_point = 0.0_dp
    real(dp) :: cdo_detach_point = 0.0_dp
    character(len=str_len) :: isin = ""
    real(dp) :: trade_price = 0.0_dp
    real(dp) :: yield = 0.0_dp
    character(len=str_len) :: payment_frequency = ""
    character(len=str_len) :: maturity_date = ""
    character(len=str_len) :: coupon_type = ""
    real(dp) :: credit_risk_weight = 0.0_dp
    character(len=str_len) :: issuer = ""
    character(len=str_len) :: fx_near_leg_fields = ""
  contains
    procedure :: configure_class
    procedure :: calc_adjusted_notional
    procedure :: calc_supervisory_duration
    procedure :: calc_maturity_factor
    procedure :: calc_supervisory_delta
    procedure :: set_time_bucket
    procedure :: is_basis_swap
    procedure :: compute_variance_units
    procedure :: calc_option_maturity
    procedure :: set_fx_dynamic
    procedure :: validate_credit_risk_weight
    procedure :: cdo_supervisory_delta
    procedure :: is_derivative
  end type trade_t

  public :: split_bond_exposure
  public :: select_derivatives

contains

  subroutine configure_class(self, class_name)
    class(trade_t), intent(inout) :: self
    character(len=*), intent(in) :: class_name
    character(len=:), allocatable :: name

    name = trim(class_name)
    self%class_name = name
    self%trade_group = ""
    self%trade_type = ""
    self%subclass = ""

    select case (name)
    case ("IRDSwap")
      self%trade_group = "IRD"
      self%trade_type = "Swap"
    case ("IRDSwaption")
      self%trade_group = "IRD"
      self%trade_type = "Option"
    case ("IRDSwapVol")
      self%trade_group = "IRD"
      self%trade_type = "Swap"
    case ("IRDFuture")
      self%trade_group = "IRD"
      self%trade_type = "Future"
    case ("Bond")
      self%trade_group = "IRD"
      self%trade_type = "Bond"
    case ("BondFuture")
      self%trade_group = "IRD"
      self%trade_type = "Future"
    case ("CDS")
      self%trade_group = "Credit"
      self%trade_type = "Single"
    case ("CDX")
      self%trade_group = "Credit"
      self%trade_type = "Index"
    case ("CDOTranche")
      self%trade_group = "Credit"
      self%trade_type = "CDO"
    case ("Commodity")
      self%trade_group = "Commodity"
    case ("CommSwap")
      self%trade_group = "Commodity"
      self%trade_type = "Swap"
    case ("CommodityForward")
      self%trade_group = "Commodity"
      self%trade_type = "Forward"
    case ("Equity")
      self%trade_group = "EQ"
    case ("EquityIndexFuture")
      self%trade_group = "EQ"
      self%trade_type = "Index Future"
    case ("EquityOptionIndex")
      self%trade_group = "EQ"
      self%trade_type = "Option"
    case ("EquityOptionSingle")
      self%trade_group = "EQ"
      self%trade_type = "Option"
    case ("FxSwap")
      self%trade_group = "FX"
      self%trade_type = "Swap"
    case ("FxForward")
      self%trade_group = "FX"
      self%trade_type = "Forward"
    case ("OtherExposure")
      self%trade_group = "OtherExposure"
      self%trade_type = "OtherExposure"
    case default
      self%trade_type = name
    end select
  end subroutine configure_class

  real(dp) function calc_adjusted_notional(self) result(value)
    class(trade_t), intent(in) :: self

    if (trim(self%trade_group) == "IRD" .or. &
        trim(self%trade_group) == "Credit") then
      value = self%notional * self%calc_supervisory_duration()
    else
      value = self%notional
    end if
  end function calc_adjusted_notional

  real(dp) function calc_supervisory_duration(self) result(value)
    class(trade_t), intent(in) :: self

    if (self%use_simplified .and. self%simplified) then
      value = self%ei - self%si
    else
      value = (exp(-0.05_dp * self%si) - exp(-0.05_dp * self%ei)) / 0.05_dp
    end if
  end function calc_supervisory_duration

  real(dp) function calc_maturity_factor(self, delivery_type) result(value)
    class(trade_t), intent(in) :: self
    character(len=*), intent(in), optional :: delivery_type
    real(dp) :: maturity_input

    maturity_input = self%ei
    if (present(delivery_type)) then
      if (trim(delivery_type) /= "Physical") maturity_input = self%si
    end if

    if (maturity_input < 1.0_dp) then
      value = sqrt(max(maturity_input, 0.0_dp))
    else
      value = 1.0_dp
    end if

    if (index(self%class_name, "Future") > 0) value = 10.0_dp / 252.0_dp
    value = max(10.0_dp / 252.0_dp, value)
  end function calc_maturity_factor

  real(dp) function calc_supervisory_delta(self, supervisory_volatility) result(value)
    class(trade_t), intent(in) :: self
    real(dp), intent(in), optional :: supervisory_volatility
    real(dp) :: lambda
    real(dp) :: option_maturity
    real(dp) :: temp
    real(dp) :: volatility
    integer :: direction

    direction = buy_sell_sign(self%buy_sell)
    if (trim(self%class_name) == "CDOTranche") then
      value = self%cdo_supervisory_delta()
      return
    end if

    if (.not. present(supervisory_volatility)) then
      value = real(direction, dp)
      return
    end if
    if (self%use_simplified .and. self%simplified) then
      value = real(direction, dp)
      return
    end if

    volatility = supervisory_volatility
    if (volatility <= 0.0_dp) error stop "calc_supervisory_delta: volatility must be positive"

    if (trim(self%trade_type) == "Option" .and. &
        index(self%class_name, "Swaption") > 0) then
      option_maturity = self%si
    else
      option_maturity = self%ei
    end if
    option_maturity = max(option_maturity, 1.0e-12_dp)

    if (self%underlying_price * self%strike_price < 0.0_dp) then
      lambda = max(0.001_dp - min(self%underlying_price, self%strike_price), 0.0_dp)
      temp = (log((self%underlying_price + lambda) / &
        (self%strike_price + lambda)) + &
        0.5_dp * volatility**2 * option_maturity) / &
        (volatility * sqrt(option_maturity))
    else
      if (self%underlying_price <= 0.0_dp .or. self%strike_price <= 0.0_dp) then
        error stop "calc_supervisory_delta: prices must be positive"
      end if
      temp = (log(self%underlying_price / self%strike_price) + &
        0.5_dp * volatility**2 * option_maturity) / &
        (volatility * sqrt(option_maturity))
    end if

    select case (uppercase(trim(self%option_type)))
    case ("CALL")
      value = real(direction, dp) * normal_cdf(temp)
    case ("PUT")
      value = -real(direction, dp) * normal_cdf(-temp)
    case default
      value = real(direction, dp)
    end select
  end function calc_supervisory_delta

  integer function set_time_bucket(self) result(bucket)
    class(trade_t), intent(in) :: self

    if (self%ei <= 1.0_dp) then
      bucket = 1
    else if (self%ei <= 5.0_dp) then
      bucket = 2
    else
      bucket = 3
    end if
  end function set_time_bucket

  logical function is_basis_swap(self) result(value)
    class(trade_t), intent(in) :: self
    character(len=:), allocatable :: leg_type

    value = .false.
    if (len_trim(self%pay_leg_type) == 0 .or. len_trim(self%rec_leg_type) == 0) return
    if (trim(self%pay_leg_type) /= trim(self%rec_leg_type)) return

    leg_type = trim(self%pay_leg_type)
    value = leg_type == "Float" .or. leg_type == "Commodity" .or. &
      leg_type == "Equity"
  end function is_basis_swap

  real(dp) function compute_variance_units(self) result(value)
    class(trade_t), intent(in) :: self

    if (abs(self%vol_strike) <= tiny(1.0_dp)) then
      error stop "compute_variance_units: vol_strike must be nonzero"
    end if
    value = self%vega_notional / (2.0_dp * self%vol_strike)
  end function compute_variance_units

  real(dp) function calc_option_maturity(self) result(value)
    class(trade_t), intent(in) :: self

    if (trim(self%delivery_type) == "Cash") then
      value = self%maturity
    else
      value = self%ei
    end if
  end function calc_option_maturity

  subroutine set_fx_dynamic(self)
    class(trade_t), intent(inout) :: self
    integer :: paying_priority
    integer :: receiving_priority

    if (trim(self%trade_group) /= "FX") return

    if (len_trim(self%ccy_pair) == 0 .or. len_trim(self%buy_sell) == 0) then
      if (len_trim(self%ccy_receiving) > 0 .and. len_trim(self%ccy_paying) > 0) then
        receiving_priority = currency_priority(self%ccy_receiving)
        paying_priority = currency_priority(self%ccy_paying)

        if (receiving_priority <= paying_priority) then
          self%ccy_pair = trim(self%ccy_receiving) // "/" // trim(self%ccy_paying)
          self%buy_sell = "Buy"
        else
          self%ccy_pair = trim(self%ccy_paying) // "/" // trim(self%ccy_receiving)
          self%buy_sell = "Sell"
        end if
      end if
    end if

    if (abs(self%notional) <= tiny(1.0_dp)) then
      if (trim(self%ccy_receiving) == trim(self%base_ccy)) then
        self%notional = self%amount_paying
      else if (trim(self%ccy_paying) == trim(self%base_ccy)) then
        self%notional = self%amount_receiving
      else
        self%notional = max(self%amount_paying, self%amount_receiving)
      end if
    end if
  end subroutine set_fx_dynamic

  logical function validate_credit_risk_weight(self) result(valid)
    class(trade_t), intent(in) :: self

    valid = self%credit_risk_weight >= 0.0_dp .and. &
      self%credit_risk_weight <= 1.0_dp
  end function validate_credit_risk_weight

  real(dp) function cdo_supervisory_delta(self) result(value)
    class(trade_t), intent(in) :: self
    real(dp) :: magnitude

    magnitude = 15.0_dp / ((1.0_dp + 14.0_dp * self%cdo_attach_point) * &
      (1.0_dp + 14.0_dp * self%cdo_detach_point))
    value = real(buy_sell_sign(self%buy_sell), dp) * magnitude
  end function cdo_supervisory_delta

  logical function is_derivative(self) result(value)
    class(trade_t), intent(in) :: self
    character(len=:), allocatable :: name
    character(len=:), allocatable :: kind

    name = uppercase(trim(self%class_name))
    kind = uppercase(trim(self%trade_type))
    value = index(name, "SWAP") > 0 .or. index(name, "FORWARD") > 0 .or. &
      index(name, "CDS") > 0 .or. index(name, "CDX") > 0 .or. &
      index(name, "CDO") > 0 .or. index(name, "FUTURE") > 0 .or. &
      index(name, "OPTION") > 0 .or. name == "OTHEREXPOSURE" .or. &
      index(kind, "SWAP") > 0 .or. index(kind, "FORWARD") > 0 .or. &
      index(kind, "FUTURE") > 0 .or. index(kind, "OPTION") > 0
  end function is_derivative

  subroutine split_bond_exposure(bond, interest_rate_trade, credit_trade)
    type(trade_t), intent(in) :: bond
    type(trade_t), intent(out) :: interest_rate_trade
    type(trade_t), intent(out) :: credit_trade
    character(len=str_len) :: interest_direction
    character(len=str_len) :: credit_direction

    if (.not. bond%validate_credit_risk_weight()) then
      error stop "split_bond_exposure: invalid credit risk weight"
    end if

    if (trim(bond%coupon_type) == "Fixed" .and. &
        uppercase(trim(bond%buy_sell)) == "BUY") then
      interest_direction = "Sell"
    else
      interest_direction = "Buy"
    end if

    if (uppercase(trim(bond%buy_sell)) == "BUY") then
      credit_direction = "Sell"
    else
      credit_direction = "Buy"
    end if

    call interest_rate_trade%configure_class("IRDSwap")
    interest_rate_trade%notional = bond%credit_risk_weight * bond%notional
    interest_rate_trade%mtm = bond%credit_risk_weight * bond%mtm
    interest_rate_trade%currency = bond%currency
    interest_rate_trade%si = bond%si
    interest_rate_trade%ei = bond%ei
    interest_rate_trade%buy_sell = interest_direction

    call credit_trade%configure_class("CDS")
    credit_trade%notional = bond%credit_risk_weight * bond%notional
    credit_trade%mtm = bond%credit_risk_weight * bond%mtm
    credit_trade%currency = bond%currency
    credit_trade%si = bond%si
    credit_trade%ei = bond%ei
    credit_trade%buy_sell = credit_direction
    credit_trade%subclass = "AA"
    credit_trade%reference_entity = bond%issuer
  end subroutine split_bond_exposure

  subroutine select_derivatives(trades, derivatives)
    type(trade_t), intent(in) :: trades(:)
    type(trade_t), allocatable, intent(out) :: derivatives(:)
    integer :: count_derivatives
    integer :: i
    integer :: output_index

    count_derivatives = 0
    do i = 1, size(trades)
      if (trades(i)%is_derivative()) count_derivatives = count_derivatives + 1
    end do

    allocate(derivatives(count_derivatives))
    output_index = 0
    do i = 1, size(trades)
      if (trades(i)%is_derivative()) then
        output_index = output_index + 1
        derivatives(output_index) = trades(i)
      end if
    end do
  end subroutine select_derivatives

  pure integer function buy_sell_sign(buy_sell) result(direction)
    character(len=*), intent(in) :: buy_sell

    if (uppercase(trim(buy_sell)) == "SELL") then
      direction = -1
    else
      direction = 1
    end if
  end function buy_sell_sign

  pure integer function currency_priority(currency) result(priority)
    character(len=*), intent(in) :: currency
    character(len=3), parameter :: currencies(9) = &
      ["EUR", "JPY", "USD", "CZK", "HRK", "HUF", "BAM", "RSD", "RUB"]
    integer :: i

    priority = 100
    do i = 1, size(currencies)
      if (uppercase(trim(currency)) == currencies(i)) then
        priority = i
        return
      end if
    end do
  end function currency_priority

  pure function uppercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: code
    integer :: i

    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        value(i:i) = achar(code - iachar('a') + iachar('A'))
      end if
    end do
  end function uppercase

end module trading_trades
