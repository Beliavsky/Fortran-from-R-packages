module saccr_core
  use trading, only : dp, str_len, trade_t, csa_t, collateral_t
  use saccr_types, only : supervisory_record_t, single_trade_addon_t, &
    replacement_cost_t, fx_hedge_result_t
  use saccr_supervisory, only : supervisory_option_volatility
  implicit none
  private

  public :: calculate_factor_multiplier
  public :: single_trade_addon
  public :: calc_ead
  public :: calc_pfe
  public :: calc_rc
  public :: determine_ccr_methodology
  public :: apply_fx_hedge

contains

  pure real(dp) function calculate_factor_multiplier(hedging_set_name) result(value)
    character(len=*), intent(in) :: hedging_set_name
    character(len=:), allocatable :: name

    name = uppercase(trim(adjustl(hedging_set_name)))
    if (index(name, "VOL_") == 1) then
      value = 5.0_dp
    else if (index(name, "BASIS_") == 1) then
      value = 0.5_dp
    else
      value = 1.0_dp
    end if
  end function calculate_factor_multiplier

  subroutine single_trade_addon(trade, records, result, maturity_factor, simplified, hedging_set)
    type(trade_t), intent(in) :: trade
    type(supervisory_record_t), intent(in) :: records(:)
    type(single_trade_addon_t), intent(out) :: result
    real(dp), intent(in), optional :: maturity_factor
    logical, intent(in), optional :: simplified
    character(len=*), intent(in), optional :: hedging_set
    type(trade_t) :: work
    character(len=str_len) :: lookup_asset
    character(len=str_len) :: lookup_subclass
    logical :: is_simplified
    logical :: found
    logical :: is_future
    real(dp) :: volatility

    work = trade
    is_simplified = .false.
    if (present(simplified)) is_simplified = simplified
    work%simplified = is_simplified
    work%use_simplified = .true.

    result = single_trade_addon_t()
    result%external_id = work%external_id
    result%asset_class = work%trade_group
    if (present(hedging_set)) result%hedging_set = hedging_set
    result%adjusted_notional = work%calc_adjusted_notional()

    is_future = index(uppercase(trim(work%class_name)), "FUTURE") > 0
    if (present(maturity_factor) .and. .not. is_future) then
      result%maturity_factor = maturity_factor
    else
      result%maturity_factor = work%calc_maturity_factor()
    end if

    if (same_text(work%trade_type, "Option")) then
      call option_lookup_key(work, lookup_asset, lookup_subclass)
      volatility = supervisory_option_volatility( &
        records, lookup_asset, lookup_subclass, found)
      if (.not. found) then
        if (same_text(work%trade_group, "Commodity")) then
          if (same_text(work%commodity_type, "Electricity")) then
            volatility = 1.5_dp
          else
            volatility = 0.7_dp
          end if
        else
          volatility = 1.5_dp
        end if
      end if
      result%volatility = volatility
      result%has_volatility = .true.
      result%supervisory_delta = work%calc_supervisory_delta(volatility)
    else
      result%supervisory_delta = work%calc_supervisory_delta()
    end if

    result%effective_notional = result%supervisory_delta * &
      result%adjusted_notional * result%maturity_factor
  end subroutine single_trade_addon

  pure real(dp) function calc_ead(rc, pfe) result(ead)
    real(dp), intent(in) :: rc
    real(dp), intent(in) :: pfe

    ead = 1.4_dp * (rc + pfe)
  end function calc_ead

  real(dp) function calc_pfe(v_c, addon_aggregate, v, simplified) result(pfe)
    real(dp), intent(in) :: v_c
    real(dp), intent(in) :: addon_aggregate
    real(dp), intent(in), optional :: v
    logical, intent(in), optional :: simplified
    real(dp) :: adjusted_v_c
    real(dp) :: multiplier
    logical :: use_simplified

    adjusted_v_c = v_c
    if (present(v)) adjusted_v_c = min(adjusted_v_c, v)

    use_simplified = .false.
    if (present(simplified)) use_simplified = simplified

    if (addon_aggregate <= 0.0_dp) then
      pfe = 0.0_dp
      return
    end if

    if (use_simplified .or. adjusted_v_c >= 0.0_dp) then
      multiplier = 1.0_dp
    else
      multiplier = min(1.0_dp, 0.05_dp + 0.95_dp * &
        exp(adjusted_v_c / (1.9_dp * addon_aggregate)))
    end if
    pfe = multiplier * addon_aggregate
  end function calc_pfe

  subroutine calc_rc(trades, result, csa, collaterals, simplified, ignore_margin)
    type(trade_t), intent(in) :: trades(:)
    type(replacement_cost_t), intent(out) :: result
    type(csa_t), intent(in), optional :: csa
    type(collateral_t), intent(in), optional :: collaterals(:)
    logical, intent(in), optional :: simplified
    logical, intent(in), optional :: ignore_margin
    real(dp) :: current_collateral
    real(dp) :: im_counterparty
    real(dp) :: mta_counterparty
    real(dp) :: scale
    real(dp) :: threshold_counterparty
    logical :: margin_ignored
    logical :: use_simplified
    integer :: i

    result = replacement_cost_t()
    if (size(trades) > 0) result%v = sum(trades%mtm)

    use_simplified = .false.
    if (present(simplified)) use_simplified = simplified
    margin_ignored = .false.
    if (present(ignore_margin)) margin_ignored = ignore_margin

    if (.not. present(csa) .or. margin_ignored) then
      result%v_c = result%v
      result%rc = max(result%v_c, 0.0_dp)
      return
    end if

    threshold_counterparty = csa%threshold_counterparty
    mta_counterparty = csa%minimum_transfer_counterparty
    im_counterparty = csa%initial_margin_counterparty

    if (same_text(csa%values_type, "Percentage") .or. &
        same_text(csa%values_type, "Perc")) then
      scale = abs(result%v)
      threshold_counterparty = threshold_counterparty * scale
      mta_counterparty = mta_counterparty * scale
      im_counterparty = im_counterparty * scale
    end if

    current_collateral = 0.0_dp
    if (present(collaterals)) then
      do i = 1, size(collaterals)
        if (.not. same_text(collaterals(i)%csa_id, csa%id)) cycle
        if (same_text(collaterals(i)%collateral_type, "ICA")) then
          current_collateral = current_collateral + collaterals(i)%amount
          im_counterparty = collaterals(i)%amount
        else if (same_text(collaterals(i)%collateral_type, "VariationMargin")) then
          current_collateral = current_collateral + collaterals(i)%amount
        end if
      end do
    end if

    if (use_simplified) then
      result%v_c = result%v
      result%rc = max(result%v_c, threshold_counterparty + mta_counterparty, 0.0_dp)
    else
      result%v_c = result%v - current_collateral
      result%rc = max(result%v_c, threshold_counterparty + &
        mta_counterparty - im_counterparty, 0.0_dp)
    end if
  end subroutine calc_rc

  function determine_ccr_methodology(trades, total_assets) result(methodology)
    type(trade_t), intent(in) :: trades(:)
    real(dp), intent(in) :: total_assets
    character(len=str_len) :: methodology
    real(dp) :: derivatives_mv
    integer :: i

    if (total_assets <= 0.0_dp) then
      error stop "determine_ccr_methodology: total_assets must be positive"
    end if

    derivatives_mv = 0.0_dp
    do i = 1, size(trades)
      if (trades(i)%is_derivative()) derivatives_mv = derivatives_mv + abs(trades(i)%mtm)
    end do

    if (derivatives_mv / total_assets < 0.05_dp .and. &
        derivatives_mv < 100000000.0_dp) then
      methodology = "Original Exposure Method"
    else if (derivatives_mv / total_assets < 0.10_dp .and. &
        derivatives_mv < 300000000.0_dp) then
      methodology = "Simplified SA-CCR"
    else
      methodology = "SA-CCR"
    end if
  end function determine_ccr_methodology

  subroutine apply_fx_hedge(ead, addon, rwa_fx_counterparty, rwa_cds_counterparty, &
      ead_cds, approach, result, protection_percentage)
    real(dp), intent(in) :: ead
    real(dp), intent(in) :: addon
    real(dp), intent(in) :: rwa_fx_counterparty
    real(dp), intent(in) :: rwa_cds_counterparty
    real(dp), intent(in) :: ead_cds
    character(len=*), intent(in) :: approach
    type(fx_hedge_result_t), intent(out) :: result
    real(dp), intent(in), optional :: protection_percentage
    real(dp) :: percentage

    result = fx_hedge_result_t()
    select case (uppercase(trim(approach)))
    case ("CURRENT")
      result%unprotected_amount = max(ead - ead_cds, 0.0_dp)
      result%protected_amount = min(max(ead_cds, 0.0_dp), ead)
    case ("TECHNICALAMENDMENT")
      if (addon <= 0.0_dp) then
        result%multiplier = 0.05_dp
      else
        result%multiplier = 0.05_dp + 0.95_dp * exp(-ead_cds / (1.9_dp * addon))
      end if
      result%unprotected_amount = result%multiplier * ead
      result%protected_amount = ead - result%unprotected_amount
    case ("CAPPEDPROTECTION")
      if (.not. present(protection_percentage)) then
        error stop "apply_fx_hedge: protection_percentage is required"
      end if
      percentage = protection_percentage
      if (percentage <= 0.0_dp .or. percentage > 1.0_dp) then
        error stop "apply_fx_hedge: protection_percentage must be in (0,1]"
      end if
      if (addon <= 0.0_dp) then
        result%multiplier = 0.05_dp
      else
        result%multiplier = 0.05_dp + 0.95_dp * &
          exp(-ead_cds / percentage / (1.9_dp * addon))
      end if
      result%unprotected_amount = (1.0_dp - percentage) * ead + &
        percentage * result%multiplier * ead
      result%protected_amount = ead - result%unprotected_amount
    case default
      error stop "apply_fx_hedge: unsupported hedging approach"
    end select

    result%fx_rwa = result%unprotected_amount * rwa_fx_counterparty
    result%cds_rwa = result%protected_amount * rwa_cds_counterparty
  end subroutine apply_fx_hedge

  subroutine option_lookup_key(trade, asset_class, subclass)
    type(trade_t), intent(in) :: trade
    character(len=str_len), intent(out) :: asset_class
    character(len=str_len), intent(out) :: subclass

    asset_class = trade%trade_group
    subclass = trade%subclass
    select case (uppercase(trim(trade%trade_group)))
    case ("IRD", "FX", "OTHEREXPOSURE")
      subclass = ""
    case ("CREDIT")
      if (same_text(trade%class_name, "CDS")) then
        asset_class = "CreditSingle"
      else
        asset_class = "CreditIndex"
      end if
    case ("EQ")
      if (.not. same_text(trade%subclass, "Index")) subclass = ""
    case ("COMMODITY")
      if (same_text(trade%commodity_type, "Electricity")) then
        subclass = "Electricity"
      else
        subclass = "Other"
      end if
    end select
  end subroutine option_lookup_key

  pure logical function same_text(left, right) result(equal)
    character(len=*), intent(in) :: left
    character(len=*), intent(in) :: right

    equal = uppercase(trim(adjustl(left))) == uppercase(trim(adjustl(right)))
  end function same_text

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

end module saccr_core
