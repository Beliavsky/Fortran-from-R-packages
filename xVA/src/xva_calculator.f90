module xva_calculator_mod
  use trading, only : dp, trade_t, csa_t, collateral_t, curve_t
  use xva_types, only : simulation_data_t, regulatory_data_t, xva_result_t, &
    exposure_profile_t, ead_result_t, cva_sensitivity_t, &
    supervisory_cva_data_t
  use xva_core, only : generate_time_grid, calc_pd, calc_va, &
    calc_effective_maturity, calc_kva, aggregate_geometric_pd
  use xva_exposure, only : calc_simulated_exposure
  use xva_regulatory, only : calc_ead_regulatory, calc_cva_capital
  use xva_supervisory, only : default_supervisory_cva_data
  use xva_math, only : normal_pdf, inverse_normal_cdf, same_text
  implicit none
  private

  public :: xva_calculator

contains

  subroutine xva_calculator(trades, csa, collateral, sim_data, reg_data, &
      credit_curve_processing_organization, credit_curve_counterparty, &
      funding_curve, spot_rates, counterparty_lgd, processing_organization_lgd, &
      no_simulations, result, seed, supervisory_data)
    type(trade_t), intent(in) :: trades(:)
    type(csa_t), intent(in) :: csa
    type(collateral_t), intent(in) :: collateral
    type(simulation_data_t), intent(in) :: sim_data
    type(regulatory_data_t), intent(in) :: reg_data
    type(curve_t), intent(in) :: credit_curve_processing_organization
    type(curve_t), intent(in) :: credit_curve_counterparty
    type(curve_t), intent(in) :: funding_curve
    type(curve_t), intent(in) :: spot_rates
    real(dp), intent(in) :: counterparty_lgd
    real(dp), intent(in) :: processing_organization_lgd
    logical, intent(in) :: no_simulations
    type(xva_result_t), intent(out) :: result
    integer, intent(in), optional :: seed
    type(supervisory_cva_data_t), intent(in), optional :: supervisory_data
    type(supervisory_cva_data_t) :: supervisory
    type(cva_sensitivity_t) :: sensitivities
    type(ead_result_t) :: ead_result
    type(exposure_profile_t) :: bumped_exposure
    type(curve_t) :: bumped_curve
    type(trade_t), allocatable :: work_trades(:)
    real(dp), allocatable :: counterparty_spread(:)
    real(dp), allocatable :: processing_organization_spread(:)
    real(dp), allocatable :: funding_spread(:)
    real(dp), allocatable :: spot_curve(:)
    real(dp), allocatable :: discount_factors(:)
    real(dp), allocatable :: pd_counterparty(:)
    real(dp), allocatable :: pd_processing_organization(:)
    real(dp), allocatable :: pd_funding(:)
    real(dp), allocatable :: bumped_spot_curve(:)
    real(dp), allocatable :: bumped_discount_factors(:)
    real(dp), allocatable :: time_points(:)
    real(dp) :: base_cva
    real(dp) :: bumped_cva
    real(dp) :: maturity
    real(dp) :: negative_exposure
    real(dp) :: pd_aggregate
    real(dp) :: positive_exposure
    real(dp), parameter :: basis_point = 0.0001_dp
    integer :: effective_seed
    integer :: i
    integer :: n
    integer :: years

    result = xva_result_t()
    if (size(trades) == 0) error stop "xva_calculator: at least one trade is required"
    if (no_simulations .and. same_text(reg_data%ccr_framework, "IMM")) then
      error stop "xva_calculator: IMM requires simulated exposures"
    end if
    if (no_simulations .and. same_text(reg_data%cva_framework, "SA-CVA")) then
      error stop "xva_calculator: SA-CVA requires simulated exposures"
    end if
    if (counterparty_lgd <= 0.0_dp .or. processing_organization_lgd <= 0.0_dp) then
      error stop "xva_calculator: LGD values must be positive"
    end if

    effective_seed = 104729
    if (present(seed)) effective_seed = seed
    if (present(supervisory_data)) then
      supervisory = supervisory_data
    else
      call default_supervisory_cva_data(supervisory)
    end if

    maturity = maxval(trades%ei)
    call generate_time_grid(csa, maturity, time_points)
    n = size(time_points)
    allocate(spot_curve(n), counterparty_spread(n), &
      processing_organization_spread(n), funding_spread(n), discount_factors(n))
    call spot_rates%interpolate(time_points, spot_curve, "linear")
    call credit_curve_counterparty%interpolate(time_points, counterparty_spread)
    call credit_curve_processing_organization%interpolate( &
      time_points, processing_organization_spread)
    call funding_curve%interpolate(time_points, funding_spread)
    discount_factors = exp(-time_points * spot_curve)

    call calc_pd(counterparty_spread, counterparty_lgd, time_points, pd_counterparty)
    call calc_pd(processing_organization_spread, processing_organization_lgd, &
      time_points, pd_processing_organization)
    call calc_pd(funding_spread, 1.0_dp, time_points, pd_funding)

    allocate(work_trades(size(trades)))
    work_trades = trades
    if (.not. no_simulations) then
      call calc_simulated_exposure(discount_factors, time_points, spot_curve, csa, &
        work_trades, sim_data, reg_data%ccr_framework, result%exposure_profile, &
        effective_seed)
      result%has_simulated_values = .true.
    else
      call allocate_zero_exposure_profile(n, result%exposure_profile)
    end if

    call calc_ead_regulatory(work_trades, reg_data%ccr_framework, &
      reg_data%sa_ccr_simplified, csa, collateral, &
      result%exposure_profile%eee, time_points, ead_result)
    result%ead = ead_result%ead_value
    result%effective_maturity = calc_effective_maturity(work_trades, time_points, &
      reg_data%ccr_framework, result%exposure_profile%ee)

    if (.not. no_simulations) then
      result%cva_simulated = calc_va(result%exposure_profile%ee, &
        discount_factors, pd_counterparty, counterparty_lgd)
      result%dva_simulated = calc_va(result%exposure_profile%nee, &
        discount_factors, pd_processing_organization, &
        processing_organization_lgd)
      result%fca_simulated = calc_va(result%exposure_profile%ee, &
        discount_factors, pd_funding)
      result%fba_simulated = calc_va(result%exposure_profile%nee, &
        discount_factors, pd_funding)
      result%mva_simulated = result%fca_simulated * 2.0_dp * &
        sqrt(reg_data%mva_days / (250.0_dp * maturity)) * &
        inverse_normal_cdf(reg_data%mva_percentile) / normal_pdf(0.0_dp)
    end if

    if (same_text(reg_data%cva_framework, "SA-CVA")) then
      base_cva = result%cva_simulated
      call calculate_credit_spread_sensitivities(credit_curve_counterparty, &
        time_points, discount_factors, result%exposure_profile%ee, &
        counterparty_lgd, pd_counterparty, base_cva, sensitivities)

      allocate(sensitivities%interest_rate_delta(size(spot_rates%rates)), &
        sensitivities%interest_rate_tenors(size(spot_rates%tenors)))
      sensitivities%interest_rate_tenors = spot_rates%tenors
      do i = 1, size(spot_rates%rates)
        bumped_curve = spot_rates
        bumped_curve%rates(i) = bumped_curve%rates(i) + basis_point
        allocate(bumped_spot_curve(n), bumped_discount_factors(n))
        call bumped_curve%interpolate(time_points, bumped_spot_curve, "linear")
        bumped_discount_factors = exp(-time_points * bumped_spot_curve)
        work_trades = trades
        call calc_simulated_exposure(bumped_discount_factors, time_points, &
          bumped_spot_curve, csa, work_trades, sim_data, reg_data%ccr_framework, &
          bumped_exposure, effective_seed)
        bumped_cva = calc_va(bumped_exposure%ee, bumped_discount_factors, &
          pd_counterparty, counterparty_lgd)
        sensitivities%interest_rate_delta(i) = -(bumped_cva - base_cva) / basis_point
        deallocate(bumped_spot_curve, bumped_discount_factors)
      end do
      call calc_cva_capital(work_trades, result%ead, reg_data, supervisory, &
        result%effective_maturity, result%cva_capital, sensitivities)
    else
      call calc_cva_capital(work_trades, result%ead, reg_data, supervisory, &
        result%effective_maturity, result%cva_capital)
    end if

    result%kva = calc_kva(result%ead, reg_data, result%effective_maturity, &
      result%cva_capital%total_charge)

    if (same_text(reg_data%ccr_framework, "SA-CCR")) then
      positive_exposure = max(ead_result%exposure%replacement_cost%v_c + &
        ead_result%exposure%addon%addon, 0.0_dp)
      if (ead_result%exposure%replacement_cost%v_c - &
          ead_result%exposure%addon%addon > 0.0_dp) then
        negative_exposure = 0.0_dp
      else if (ead_result%exposure%replacement_cost%v_c < 0.0_dp) then
        negative_exposure = ead_result%exposure%replacement_cost%v_c - &
          ead_result%exposure%addon%addon
      else
        negative_exposure = -ead_result%exposure%replacement_cost%v_c - &
          ead_result%exposure%addon%addon
      end if

      years = ceiling(result%effective_maturity)
      pd_aggregate = aggregate_geometric_pd(reg_data%pd_counterparty, years)
      result%cva_saccr = -positive_exposure * pd_aggregate * counterparty_lgd
      pd_aggregate = aggregate_geometric_pd( &
        reg_data%pd_processing_organization, years)
      result%dva_saccr = -negative_exposure * pd_aggregate * &
        processing_organization_lgd
      pd_aggregate = aggregate_geometric_pd(reg_data%pd_counterparty, years, &
        reg_data%pd_funding)
      result%fca_saccr = -positive_exposure * pd_aggregate
      result%fba_saccr = -negative_exposure * pd_aggregate
      result%mva_saccr = result%fca_saccr * 2.0_dp * &
        sqrt(reg_data%mva_days / (250.0_dp * maturity)) * &
        inverse_normal_cdf(reg_data%mva_percentile) / normal_pdf(0.0_dp)
      result%has_saccr_values = .true.
    end if
  end subroutine xva_calculator

  subroutine calculate_credit_spread_sensitivities(curve, time_points, &
      discount_factors, exposure, lgd, base_pd, base_cva, sensitivities)
    type(curve_t), intent(in) :: curve
    real(dp), intent(in) :: time_points(:)
    real(dp), intent(in) :: discount_factors(:)
    real(dp), intent(in) :: exposure(:)
    real(dp), intent(in) :: lgd
    real(dp), intent(in) :: base_pd(:)
    real(dp), intent(in) :: base_cva
    type(cva_sensitivity_t), intent(inout) :: sensitivities
    type(curve_t) :: bumped_curve
    real(dp), allocatable :: bumped_spread(:)
    real(dp), allocatable :: bumped_pd(:)
    real(dp) :: bumped_cva
    real(dp), parameter :: basis_point = 0.0001_dp
    integer :: i

    allocate(sensitivities%credit_spread_delta(size(curve%rates)), &
      sensitivities%credit_spread_tenors(size(curve%tenors)), &
      bumped_spread(size(time_points)))
    sensitivities%credit_spread_tenors = curve%tenors
    do i = 1, size(curve%rates)
      bumped_curve = curve
      bumped_curve%rates(i) = bumped_curve%rates(i) + 1.0_dp
      call bumped_curve%interpolate(time_points, bumped_spread)
      call calc_pd(bumped_spread, lgd, time_points, bumped_pd)
      bumped_pd = max(bumped_pd, base_pd)
      bumped_cva = calc_va(exposure, discount_factors, bumped_pd, lgd)
      sensitivities%credit_spread_delta(i) = &
        -(bumped_cva - base_cva) / basis_point
    end do
  end subroutine calculate_credit_spread_sensitivities

  subroutine allocate_zero_exposure_profile(n, profile)
    integer, intent(in) :: n
    type(exposure_profile_t), intent(out) :: profile

    allocate(profile%ee_uncollateralized(n), profile%nee_uncollateralized(n), &
      profile%pfe_uncollateralized(n), profile%ee(n), profile%nee(n), &
      profile%pfe(n), profile%eee(n))
    profile%ee_uncollateralized = 0.0_dp
    profile%nee_uncollateralized = 0.0_dp
    profile%pfe_uncollateralized = 0.0_dp
    profile%ee = 0.0_dp
    profile%nee = 0.0_dp
    profile%pfe = 0.0_dp
    profile%eee = 0.0_dp
  end subroutine allocate_zero_exposure_profile

end module xva_calculator_mod
