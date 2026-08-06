module xva_regulatory
  use trading, only : dp, trade_t, csa_t, collateral_t
  use saccr, only : exposure_result_t, calculate_exposure
  use xva_types, only : ead_result_t, regulatory_data_t, &
    supervisory_cva_data_t, cva_sensitivity_t, cva_capital_result_t
  use xva_core, only : calc_ngr, cem_addon_factor, is_eligible_currency, &
    is_investment_grade
  use xva_supervisory, only : rating_weight, sector_risk_weight
  use xva_math, only : inverse_normal_cdf, quadratic_form, same_text
  implicit none
  private

  public :: calc_ead_regulatory
  public :: calc_cva_capital

contains

  subroutine calc_ead_regulatory(trades, framework, sa_ccr_simplified, csa, &
      collateral, eee, time_points, result)
    type(trade_t), intent(in) :: trades(:)
    character(len=*), intent(in) :: framework
    character(len=*), intent(in) :: sa_ccr_simplified
    type(csa_t), intent(in) :: csa
    type(collateral_t), intent(in) :: collateral
    real(dp), intent(in) :: eee(:)
    real(dp), intent(in) :: time_points(:)
    type(ead_result_t), intent(out) :: result
    type(collateral_t) :: collaterals(1)
    real(dp), allocatable :: mtm_values(:)
    real(dp) :: addon_amount
    real(dp) :: initial_margin
    real(dp) :: ngr
    integer :: i
    integer :: points_within_year
    logical :: oem
    logical :: simplified

    result = ead_result_t()
    if (size(trades) == 0) return

    if (same_text(framework, "CEM")) then
      allocate(mtm_values(size(trades)))
      mtm_values = trades%mtm
      ngr = calc_ngr(mtm_values)
      addon_amount = 0.0_dp
      do i = 1, size(trades)
        addon_amount = addon_amount + abs(trades(i)%notional) * &
          cem_addon_factor(trades(i)%ei)
      end do
      addon_amount = (0.4_dp + 0.6_dp * ngr) * addon_amount
      initial_margin = max(csa%initial_margin_counterparty, 0.0_dp)
      result%ead_value = max(max(sum(mtm_values), 0.0_dp) + &
        addon_amount - initial_margin, 0.0_dp)

    else if (same_text(framework, "SA-CCR")) then
      simplified = same_text(sa_ccr_simplified, "simplified") .or. &
        same_text(sa_ccr_simplified, "OEM")
      oem = same_text(sa_ccr_simplified, "OEM")
      collaterals(1) = collateral
      call calculate_exposure(trades, result%exposure, csa=csa, &
        collaterals=collaterals, simplified=simplified, oem=oem)
      result%ead_value = result%exposure%ead
      result%has_exposure = .true.

    else if (same_text(framework, "IMM")) then
      if (size(eee) /= size(time_points)) then
        error stop "calc_ead_regulatory: EEE and time_points sizes differ"
      end if
      points_within_year = count(time_points <= 1.0_dp)
      if (points_within_year == 0) then
        error stop "calc_ead_regulatory: IMM requires time points within one year"
      end if
      result%ead_value = 1.4_dp * sum(eee, mask=time_points <= 1.0_dp) / &
        real(points_within_year, dp)
    else
      error stop "calc_ead_regulatory: framework must be CEM, SA-CCR, or IMM"
    end if
  end subroutine calc_ead_regulatory

  subroutine calc_cva_capital(trades, ead, reg_data, supervisory, &
      effective_maturity, result, sensitivities)
    type(trade_t), intent(in) :: trades(:)
    real(dp), intent(in) :: ead
    type(regulatory_data_t), intent(in) :: reg_data
    type(supervisory_cva_data_t), intent(in) :: supervisory
    real(dp), intent(in) :: effective_maturity
    type(cva_capital_result_t), intent(out) :: result
    type(cva_sensitivity_t), intent(in), optional :: sensitivities
    real(dp), allocatable :: mapped(:)
    real(dp), allocatable :: weighted(:)
    real(dp) :: discount_factor
    real(dp) :: regulatory_weight
    real(dp) :: supervisory_weight
    logical :: found

    result = cva_capital_result_t()
    regulatory_weight = rating_weight(supervisory, &
      reg_data%counterparty_rating, found)
    if (.not. found) error stop "calc_cva_capital: unknown counterparty rating"
    supervisory_weight = sector_risk_weight(supervisory, &
      reg_data%counterparty_sector, &
      is_investment_grade(reg_data%counterparty_rating), found)
    if (.not. found) error stop "calc_cva_capital: unknown counterparty sector"

    if (same_text(reg_data%cva_framework, "STD-CVA")) then
      if (effective_maturity <= 0.0_dp) then
        error stop "calc_cva_capital: effective maturity must be positive"
      end if
      discount_factor = (1.0_dp - exp(-0.05_dp * effective_maturity)) / &
        (effective_maturity * 0.05_dp)
      result%total_charge = inverse_normal_cdf(0.99_dp) * abs( &
        regulatory_weight * ead * discount_factor * effective_maturity)

    else if (same_text(reg_data%cva_framework, "BA-CVA")) then
      result%total_charge = ead / 1.4_dp * supervisory_weight * &
        effective_maturity

    else if (same_text(reg_data%cva_framework, "SA-CVA")) then
      if (.not. present(sensitivities)) then
        error stop "calc_cva_capital: sensitivities are required for SA-CVA"
      end if
      call validate_sensitivities(sensitivities)
      if (size(trades) == 0) error stop "calc_cva_capital: trades are required"

      if (is_eligible_currency(trades(1)%currency)) then
        call map_sensitivities(sensitivities%interest_rate_delta, &
          sensitivities%interest_rate_tenors, supervisory%ir_eligible_tenors, mapped)
        allocate(weighted(size(mapped)))
        weighted = mapped * supervisory%ir_risk_weight_eligible
        result%interest_rate_charge = reg_data%sa_cva_multiplier * &
          sqrt(max(quadratic_form(weighted, &
          supervisory%ir_eligible_correlation), 0.0_dp))
      else
        call map_sensitivities(sensitivities%interest_rate_delta, &
          sensitivities%interest_rate_tenors, supervisory%ir_other_tenors, mapped)
        allocate(weighted(size(mapped)))
        weighted = mapped * supervisory%ir_risk_weight_other
        result%interest_rate_charge = reg_data%sa_cva_multiplier * &
          sqrt(max(quadratic_form(weighted, &
          supervisory%ir_other_correlation), 0.0_dp))
      end if

      deallocate(mapped, weighted)
      call map_sensitivities(sensitivities%credit_spread_delta, &
        sensitivities%credit_spread_tenors, supervisory%cs_tenors, mapped)
      result%credit_spread_charge = reg_data%sa_cva_multiplier * &
        supervisory_weight * sqrt(max(quadratic_form(mapped, &
        supervisory%cs_tenor_correlation), 0.0_dp))
      result%total_charge = result%interest_rate_charge + &
        result%credit_spread_charge
    else
      error stop "calc_cva_capital: unsupported CVA framework"
    end if
  end subroutine calc_cva_capital

  subroutine validate_sensitivities(sensitivities)
    type(cva_sensitivity_t), intent(in) :: sensitivities

    if (.not. allocated(sensitivities%interest_rate_delta) .or. &
        .not. allocated(sensitivities%interest_rate_tenors) .or. &
        .not. allocated(sensitivities%credit_spread_delta) .or. &
        .not. allocated(sensitivities%credit_spread_tenors)) then
      error stop "validate_sensitivities: sensitivity arrays are not allocated"
    end if
    if (size(sensitivities%interest_rate_delta) /= &
        size(sensitivities%interest_rate_tenors)) then
      error stop "validate_sensitivities: IR sensitivity sizes differ"
    end if
    if (size(sensitivities%credit_spread_delta) /= &
        size(sensitivities%credit_spread_tenors)) then
      error stop "validate_sensitivities: CS sensitivity sizes differ"
    end if
  end subroutine validate_sensitivities

  subroutine map_sensitivities(values, tenors, supervisory_tenors, mapped)
    real(dp), intent(in) :: values(:)
    real(dp), intent(in) :: tenors(:)
    real(dp), intent(in) :: supervisory_tenors(:)
    real(dp), allocatable, intent(out) :: mapped(:)
    integer :: bucket
    integer :: i
    integer :: j

    if (size(values) /= size(tenors)) then
      error stop "map_sensitivities: values and tenors sizes differ"
    end if
    allocate(mapped(size(supervisory_tenors)))
    mapped = 0.0_dp
    do i = 1, size(values)
      bucket = 0
      do j = 1, size(supervisory_tenors)
        if (supervisory_tenors(j) <= tenors(i)) bucket = j
      end do
      if (bucket > 0) mapped(bucket) = mapped(bucket) + values(i)
    end do
  end subroutine map_sensitivities

end module xva_regulatory
