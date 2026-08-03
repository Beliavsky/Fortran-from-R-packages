! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
module creditr
  use creditr_kinds, only: dp, creditr_ok, creditr_invalid_input, quiet_nan
  use creditr_dates, only: date_t, cds_dates_t, conventions_t, add_dates, add_conventions, make_date
  use creditr_curve, only: rate_quote_t, zero_curve_t, build_zero_curve, read_rate_quotes_csv
  use creditr_cds, only: spread_to_upfront_curve, upfront_to_spread_curve, price_cds_points
  implicit none
  private

  type, public :: cds_contract_t
    character(len=80) :: name = 'NA'
    character(len=16) :: contract = 'SNAC'
    character(len=16) :: red = 'NA'
    type(date_t) :: trade_date
    type(date_t) :: maturity
    integer :: tenor_years = 5
    real(kind=dp) :: spread_bps = 0.0_dp
    real(kind=dp) :: coupon_bps = 100.0_dp
    real(kind=dp) :: recovery = 0.4_dp
    character(len=3) :: currency = 'USD'
    real(kind=dp) :: notional = 1.0e7_dp
    logical :: use_maturity = .false.
  end type cds_contract_t

  type, public :: cds_result_t
    type(cds_contract_t) :: contract
    type(cds_dates_t) :: dates
    type(conventions_t) :: conventions
    real(kind=dp) :: principal = 0.0_dp
    real(kind=dp) :: upfront = 0.0_dp
    real(kind=dp) :: accrual = 0.0_dp
    real(kind=dp) :: probability_default = 0.0_dp
    real(kind=dp) :: price = 100.0_dp
    real(kind=dp) :: spread_dv01 = 0.0_dp
    real(kind=dp) :: ir_dv01 = 0.0_dp
    real(kind=dp) :: recovery_risk_01 = 0.0_dp
    real(kind=dp) :: cs10 = 0.0_dp
    real(kind=dp) :: hazard_rate = 0.0_dp
    integer :: status = creditr_ok
  end type cds_result_t

  public :: dp, creditr_ok, creditr_invalid_input
  public :: date_t, cds_dates_t, conventions_t, rate_quote_t, zero_curve_t
  public :: make_date, add_dates, add_conventions, build_zero_curve, read_rate_quotes_csv
  public :: price_cds, cds, spread_to_upfront, upfront_to_spread
  public :: spread_dv01, ir_dv01, rec_risk_01, cs10
  public :: spread_to_pd, pd_to_spread, implied_rr, pv01

contains

  function price_cds(contract, curve, status, quotes) result(result)
    type(cds_contract_t), intent(in) :: contract
    type(zero_curve_t), intent(in) :: curve
    type(rate_quote_t), intent(in), optional :: quotes(:)
    integer, intent(out), optional :: status
    type(cds_result_t) :: result
    type(zero_curve_t) :: bumped_curve
    type(rate_quote_t), allocatable :: bumped_quotes(:)
    real(kind=dp) :: points, bumped, h
    integer :: local_status

    result%contract = contract
    if (contract%use_maturity) then
      call add_dates(contract%trade_date, contract%currency, result%dates, maturity=contract%maturity, &
        status=local_status)
    else
      call add_dates(contract%trade_date, contract%currency, result%dates, tenor_years=contract%tenor_years, &
        status=local_status)
    end if
    call add_conventions(contract%currency, result%conventions, local_status)
    if (local_status /= creditr_ok .or. contract%notional <= 0.0_dp) then
      result%status = creditr_invalid_input
      if (present(status)) status = result%status
      return
    end if

    call price_cds_points(curve, result%dates, contract%spread_bps, contract%coupon_bps, contract%recovery, &
      .true., points, h, local_status)
    result%principal = contract%notional * points
    result%hazard_rate = h
    call price_cds_points(curve, result%dates, contract%spread_bps, contract%coupon_bps, contract%recovery, &
      .false., points, status=local_status)
    result%upfront = contract%notional * points
    result%accrual = result%upfront - result%principal
    result%price = (1.0_dp - result%principal / contract%notional) * 100.0_dp
    result%probability_default = spread_to_pd(contract%spread_bps, contract%recovery, &
      real(result%dates%end_date%serial() - contract%trade_date%serial(), dp) / 360.0_dp)

    call price_cds_points(curve, result%dates, contract%spread_bps + 1.0_dp, contract%coupon_bps, &
      contract%recovery, .false., bumped, status=local_status)
    result%spread_dv01 = contract%notional * bumped - result%upfront
    call price_cds_points(curve, result%dates, 1.1_dp * contract%spread_bps, contract%coupon_bps, &
      contract%recovery, .false., bumped, status=local_status)
    result%cs10 = contract%notional * bumped - result%upfront
    call price_cds_points(curve, result%dates, contract%spread_bps, contract%coupon_bps, &
      contract%recovery + 0.01_dp, .false., bumped, status=local_status)
    result%recovery_risk_01 = contract%notional * bumped - result%upfront

    if (present(quotes)) then
      allocate(bumped_quotes(size(quotes)))
      bumped_quotes = quotes
      bumped_quotes%rate = bumped_quotes%rate + 1.0e-4_dp
      call build_zero_curve(curve%base_date, bumped_quotes, result%conventions, bumped_curve, local_status)
    else
      bumped_curve = bump_zero_curve(curve, 1.0e-4_dp)
    end if
    call price_cds_points(bumped_curve, result%dates, contract%spread_bps, contract%coupon_bps, &
      contract%recovery, .false., bumped, status=local_status)
    result%ir_dv01 = contract%notional * bumped - result%upfront

    result%status = local_status
    if (present(status)) status = result%status
  end function price_cds

  function cds(contract, curve, status, quotes) result(result)
    type(cds_contract_t), intent(in) :: contract
    type(zero_curve_t), intent(in) :: curve
    type(rate_quote_t), intent(in), optional :: quotes(:)
    integer, intent(out), optional :: status
    type(cds_result_t) :: result
    result = price_cds(contract, curve, status, quotes)
  end function cds

  real(kind=dp) function spread_to_upfront(curve, dates, spread_bps, coupon_bps, recovery, notional, &
      is_price_clean, status) result(value)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: spread_bps, coupon_bps, recovery, notional
    logical, intent(in), optional :: is_price_clean
    integer, intent(out), optional :: status
    value = spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery, notional, &
      is_price_clean, status)
  end function spread_to_upfront

  real(kind=dp) function upfront_to_spread(curve, dates, upfront, coupon_bps, recovery, notional, &
      is_price_clean, status) result(value)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: upfront, coupon_bps, recovery, notional
    logical, intent(in), optional :: is_price_clean
    integer, intent(out), optional :: status
    logical :: clean
    clean = .false.
    if (present(is_price_clean)) clean = is_price_clean
    if (notional <= 0.0_dp) then
      value = quiet_nan()
      if (present(status)) status = creditr_invalid_input
      return
    end if
    value = upfront_to_spread_curve(curve, dates, upfront / notional, coupon_bps, recovery, clean, status)
  end function upfront_to_spread

  real(kind=dp) function spread_dv01(curve, dates, spread_bps, coupon_bps, recovery, notional) result(value)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: spread_bps, coupon_bps, recovery, notional
    value = spread_to_upfront_curve(curve, dates, spread_bps + 1.0_dp, coupon_bps, recovery, notional) - &
      spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery, notional)
  end function spread_dv01

  real(kind=dp) function cs10(curve, dates, spread_bps, coupon_bps, recovery, notional) result(value)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: spread_bps, coupon_bps, recovery, notional
    value = spread_to_upfront_curve(curve, dates, 1.1_dp * spread_bps, coupon_bps, recovery, notional) - &
      spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery, notional)
  end function cs10

  real(kind=dp) function rec_risk_01(curve, dates, spread_bps, coupon_bps, recovery, notional) result(value)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: spread_bps, coupon_bps, recovery, notional
    value = spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery + 0.01_dp, notional) - &
      spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery, notional)
  end function rec_risk_01

  real(kind=dp) function ir_dv01(curve, dates, spread_bps, coupon_bps, recovery, notional) result(value)
    type(zero_curve_t), intent(in) :: curve
    type(cds_dates_t), intent(in) :: dates
    real(kind=dp), intent(in) :: spread_bps, coupon_bps, recovery, notional
    type(zero_curve_t) :: bumped
    bumped = bump_zero_curve(curve, 1.0e-4_dp)
    value = spread_to_upfront_curve(bumped, dates, spread_bps, coupon_bps, recovery, notional) - &
      spread_to_upfront_curve(curve, dates, spread_bps, coupon_bps, recovery, notional)
  end function ir_dv01

  pure real(kind=dp) function spread_to_pd(spread_bps, recovery, time_years) result(pd)
    real(kind=dp), intent(in) :: spread_bps, recovery, time_years
    if (recovery >= 1.0_dp .or. time_years < 0.0_dp) then
      pd = quiet_nan()
    else
      pd = 1.0_dp - exp(-spread_bps / 1.0e4_dp * time_years / (1.0_dp - recovery))
    end if
  end function spread_to_pd

  pure real(kind=dp) function pd_to_spread(pd, recovery, time_years) result(spread_bps)
    real(kind=dp), intent(in) :: pd, recovery, time_years
    if (pd < 0.0_dp .or. pd >= 1.0_dp .or. time_years <= 0.0_dp) then
      spread_bps = quiet_nan()
    else
      spread_bps = 1.0e4_dp * (recovery - 1.0_dp) * log(1.0_dp - pd) / time_years
    end if
  end function pd_to_spread

  pure real(kind=dp) function implied_rr(pd, spread_bps, tenor_years) result(recovery_percent)
    real(kind=dp), intent(in) :: pd, spread_bps, tenor_years
    if (pd <= 0.0_dp .or. pd >= 1.0_dp) then
      recovery_percent = quiet_nan()
    else
      recovery_percent = 100.0_dp + (spread_bps * tenor_years / 100.0_dp) / log(1.0_dp - pd)
    end if
  end function implied_rr

  pure real(kind=dp) function pv01(principal, notional, spread_bps, coupon_bps) result(value)
    real(kind=dp), intent(in) :: principal, notional, spread_bps, coupon_bps
    if (abs(notional) <= tiny(1.0_dp) .or. abs(spread_bps - coupon_bps) <= tiny(1.0_dp)) then
      value = quiet_nan()
    else
      value = abs(principal) / notional * 1.0e4_dp / abs(spread_bps - coupon_bps)
    end if
  end function pv01

  function bump_zero_curve(curve, bump) result(out)
    type(zero_curve_t), intent(in) :: curve
    real(kind=dp), intent(in) :: bump
    type(zero_curve_t) :: out
    integer :: i
    out%base_date = curve%base_date
    allocate(out%node_serial(size(curve%node_serial)), out%log_discount(size(curve%log_discount)))
    out%node_serial = curve%node_serial
    do i = 1, size(curve%node_serial)
      out%log_discount(i) = curve%log_discount(i) - bump * &
        real(curve%node_serial(i) - curve%base_date%serial(), dp) / 365.0_dp
    end do
  end function bump_zero_curve

end module creditr
