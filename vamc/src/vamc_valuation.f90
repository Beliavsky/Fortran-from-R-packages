module vamc_valuation
  use vamc_kinds, only: dp
  use vamc_status, only: status_type, vamc_invalid_argument, vamc_dimension_error
  use vamc_dates, only: date_type, add_months, months_between, days_between, &
                        operator(==), operator(<), operator(<=), operator(>=), operator(>)
  use vamc_policy, only: policy_type, portfolio_type
  use vamc_mortality, only: mortality_table_type, mortality_factors_type, calc_mort_factors
  implicit none
  private

  type, public :: valuation_options_type
    logical :: source_compatible_timing = .true.
    logical :: source_compatible_ab_renewal = .true.
    logical :: source_compatible_maturity_vector = .true.
    logical :: source_compatible_aging_index = .true.
  end type valuation_options_type

  type, public :: projection_result_type
    real(dp) :: death_benefit = 0.0_dp
    real(dp) :: living_benefit = 0.0_dp
    real(dp) :: risk_charge = 0.0_dp
    type(policy_type) :: out_policy
  end type projection_result_type

  type, public :: policy_valuation_type
    real(dp) :: policy_value = 0.0_dp
    real(dp) :: risk_charge = 0.0_dp
  end type policy_valuation_type

  type, public :: portfolio_valuation_type
    real(dp) :: portfolio_value = 0.0_dp
    real(dp) :: portfolio_risk_charge = 0.0_dp
    real(dp), allocatable :: policy_values(:)
    real(dp), allocatable :: policy_risk_charges(:)
  end type portfolio_valuation_type

  public :: project_policy, valuate_one_policy, valuate_portfolio, age_one_policy, age_portfolio
  interface valuate_one_policy
    module procedure valuate_one_policy_2d
    module procedure valuate_one_policy_3d
  end interface valuate_one_policy
contains
  subroutine project_policy(input_policy, fund_scenario, dt, pq, survival, discount, result, options, status)
    type(policy_type), intent(in) :: input_policy
    real(dp), intent(in) :: fund_scenario(:,:), dt, pq(:), survival(:), discount(:)
    type(projection_result_type), intent(out) :: result
    type(valuation_options_type), intent(in), optional :: options
    type(status_type), intent(inout), optional :: status
    type(policy_type) :: policy
    type(valuation_options_type) :: opt
    real(dp), allocatable :: death(:), living(:), charges(:), account_history(:), partial(:)
    real(dp) :: account, raw, fee, withdrawal_guarantee, withdrawal, ag, at
    integer :: nstep, nfund, step, original_term_months, source_months
    logical :: anniversary, maturity_event, has_death, has_ab, has_ib, has_mb, has_wb
    character(len=4) :: product
    if (present(status)) call status%clear()
    opt = valuation_options_type()
    if (present(options)) opt = options
    nstep = size(fund_scenario,1)
    nfund = size(fund_scenario,2)
    if (nstep < 1 .or. nfund < 1 .or. size(pq) < nstep .or. size(survival) < nstep .or. size(discount) < nstep) then
      if (present(status)) call status%set(vamc_dimension_error, 'Scenario, mortality, and discount dimensions do not conform.')
      return
    end if
    if (.not. allocated(input_policy%fund_values) .or. .not. allocated(input_policy%fund_fees)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Policy fund arrays are not allocated.')
      return
    end if
    if (size(input_policy%fund_values) /= nfund .or. size(input_policy%fund_fees) /= nfund) then
      if (present(status)) call status%set(vamc_dimension_error, 'Policy funds do not align with fund scenario.')
      return
    end if
    policy = input_policy
    product = policy%product_type
    has_death = product(1:2) == 'DB'
    has_ab = product(1:2) == 'AB' .or. product == 'DBAB'
    has_ib = product(1:2) == 'IB' .or. product == 'DBIB'
    has_mb = product(1:2) == 'MB' .or. product == 'DBMB'
    has_wb = product(1:2) == 'WB' .or. product == 'DBWB'
    if (.not. (has_death .or. has_ab .or. has_ib .or. has_mb .or. has_wb)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Unknown variable-annuity product type.')
      return
    end if
    allocate(death(nstep), living(nstep), charges(nstep), account_history(nstep), partial(nfund))
    death = 0.0_dp
    living = 0.0_dp
    charges = 0.0_dp
    account_history = 0.0_dp
    withdrawal = 0.0_dp
    original_term_months = max(1, months_between(policy%issue_date, policy%maturity_date))
    call annuity_factors(survival, discount, nstep, ag, at)
    do step = 1, nstep
      partial = policy%fund_values * fund_scenario(step,:) * (1.0_dp - policy%fund_fees * dt)
      anniversary = policy%current_date%month == policy%issue_date%month
      fee = 0.0_dp
      if (anniversary) then
        raw = sum(partial)
        fee = raw * policy%rider_fee
        partial = partial - partial * policy%rider_fee - partial * policy%base_fee
        account = sum(partial)
        policy%fund_values = partial
        call update_guarantee(policy, product, account)
        if (has_wb) then
          withdrawal_guarantee = policy%guarantee_amount * policy%withdrawal_rate
          withdrawal = min(withdrawal_guarantee, policy%gmwb_balance)
          policy%gmwb_balance = policy%gmwb_balance - withdrawal
          policy%cumulative_withdrawal = policy%cumulative_withdrawal + withdrawal
          if (policy%gmwb_balance < 1.0e-4_dp) policy%guarantee_amount = 0.0_dp
          account = max(0.0_dp, account - withdrawal)
          if (account > 1.0e-5_dp) then
            policy%fund_values = policy%fund_values * (account / (account + withdrawal))
          else
            policy%fund_values = 0.0_dp
          end if
        end if
      else
        account = sum(partial)
        policy%fund_values = partial
      end if
      account_history(step) = account
      if (has_death) death(step) = max(0.0_dp, policy%guarantee_amount - account)
      charges(step) = fee
      if (opt%source_compatible_timing) then
        maturity_event = policy%current_date == policy%maturity_date
      else
        maturity_event = add_months(policy%current_date,1) >= policy%maturity_date .and. &
                         policy%current_date < policy%maturity_date
      end if
      if (has_ab .and. maturity_event) then
        if (policy%guarantee_amount < account) then
          policy%guarantee_amount = account
        else
          living(step) = max(0.0_dp, policy%guarantee_amount - account)
          if (account > 1.0e-5_dp) then
            policy%fund_values = policy%fund_values * policy%guarantee_amount / account
          else
            policy%fund_values = policy%guarantee_amount / real(nfund,dp)
          end if
        end if
      else if (has_ib .and. maturity_event) then
        if (ag > tiny(1.0_dp)) living(step) = max(policy%guarantee_amount * at / ag - account, 0.0_dp)
      else if (has_mb .and. maturity_event) then
        if (opt%source_compatible_maturity_vector) then
          living(step) = max(0.0_dp, maxval(policy%guarantee_amount - account_history))
        else
          living(step) = max(policy%guarantee_amount - account, 0.0_dp)
        end if
      else if (has_wb) then
        living(step) = max(0.0_dp, withdrawal - account)
      end if
      if (has_ab) then
        if (opt%source_compatible_ab_renewal) then
          source_months = days_between(policy%issue_date, policy%current_date)
          policy%maturity_date = add_months(policy%current_date, source_months)
          policy%issue_date = policy%current_date
        else if (maturity_event) then
          policy%issue_date = policy%current_date
          policy%maturity_date = add_months(policy%current_date, original_term_months)
        end if
      end if
      policy%current_date = add_months(policy%current_date, 1)
    end do
    result%death_benefit = sum(pq(1:nstep) * discount(1:nstep) * death)
    result%living_benefit = sum(survival(1:nstep) * discount(1:nstep) * living)
    result%risk_charge = sum(survival(1:nstep) * discount(1:nstep) * charges)
    result%out_policy = policy
  contains
    subroutine update_guarantee(p, product_name, current_account)
      type(policy_type), intent(inout) :: p
      character(len=4), intent(in) :: product_name
      real(dp), intent(in) :: current_account
      if (product_name == 'DBAB' .or. product_name == 'DBIB' .or. product_name == 'DBMB' .or. product_name == 'DBWB') then
        p%guarantee_amount = max(current_account, p%guarantee_amount)
      else if (product_name(3:4) == 'RU') then
        p%guarantee_amount = p%guarantee_amount * (1.0_dp + p%roll_up_rate)
      else if (product_name(3:4) == 'SU') then
        p%guarantee_amount = max(current_account, p%guarantee_amount)
      end if
    end subroutine update_guarantee

    subroutine annuity_factors(pvec, dvec, steps, guaranteed, market)
      real(dp), intent(in) :: pvec(:), dvec(:)
      integer, intent(in) :: steps
      real(dp), intent(out) :: guaranteed, market
      integer :: ny, id
      real(dp) :: pp, forward_rate
      guaranteed = 0.0_dp
      market = 0.0_dp
      pp = 1.0_dp
      ny = 0
      do while (pp > 1.0e-5_dp .and. ny < 1000)
        id = ny * 12 + 1
        if (id <= size(pvec)) then
          pp = pvec(id)
        else
          pp = 0.0_dp
        end if
        guaranteed = guaranteed + pp * exp(-0.05_dp * real(ny,dp))
        if (id < steps) then
          forward_rate = dvec(id)
        else
          forward_rate = dvec(max(1,steps-1))
        end if
        market = market + pp * exp(-forward_rate * real(ny,dp))
        ny = ny + 1
      end do
    end subroutine annuity_factors
  end subroutine project_policy

  subroutine valuate_one_policy_2d(policy, mortality_table, fund_scenario, dt, discount, valuation, options, status)
    type(policy_type), intent(in) :: policy
    type(mortality_table_type), intent(in) :: mortality_table
    real(dp), intent(in) :: fund_scenario(:,:), dt, discount(:)
    type(policy_valuation_type), intent(out) :: valuation
    type(valuation_options_type), intent(in), optional :: options
    type(status_type), intent(inout), optional :: status
    type(mortality_factors_type) :: mortality
    type(projection_result_type) :: projection
    type(status_type) :: local_status
    integer :: nstep
    if (present(status)) call status%clear()
    valuation = policy_valuation_type()
    if (policy%maturity_date <= policy%current_date) return
    nstep = months_between(policy%current_date, policy%maturity_date)
    if (nstep < 1 .or. size(fund_scenario,1) < nstep .or. size(discount) < nstep) then
      if (present(status)) call status%set(vamc_dimension_error, 'Insufficient scenario or discount horizon.')
      return
    end if
    call calc_mort_factors(policy, mortality_table, dt, mortality, local_status)
    if (.not. local_status%ok()) then
      if (present(status)) call status%set(local_status%code, local_status%message)
      return
    end if
    if (size(mortality%pq) < nstep) then
      if (present(status)) call status%set(vamc_dimension_error, 'Mortality horizon is shorter than policy horizon.')
      return
    end if
    call project_policy(policy, fund_scenario(1:nstep,:), dt, mortality%pq, mortality%p, discount, &
                        projection, options, local_status)
    if (.not. local_status%ok()) then
      if (present(status)) call status%set(local_status%code, local_status%message)
      return
    end if
    valuation%policy_value = projection%death_benefit + projection%living_benefit
    valuation%risk_charge = projection%risk_charge
  end subroutine valuate_one_policy_2d

  subroutine valuate_one_policy_3d(policy, mortality_table, fund_scenarios, dt, discount, valuation, options, status)
    type(policy_type), intent(in) :: policy
    type(mortality_table_type), intent(in) :: mortality_table
    real(dp), intent(in) :: fund_scenarios(:,:,:), dt, discount(:)
    type(policy_valuation_type), intent(out) :: valuation
    type(valuation_options_type), intent(in), optional :: options
    type(status_type), intent(inout), optional :: status
    type(policy_valuation_type) :: one
    type(status_type) :: local_status
    integer :: i
    if (present(status)) call status%clear()
    valuation = policy_valuation_type()
    if (size(fund_scenarios,1) < 1) then
      if (present(status)) call status%set(vamc_dimension_error, 'At least one scenario is required.')
      return
    end if
    do i = 1, size(fund_scenarios,1)
      call valuate_one_policy_2d(policy, mortality_table, fund_scenarios(i,:,:), dt, discount, one, options, local_status)
      if (.not. local_status%ok()) then
        if (present(status)) call status%set(local_status%code, local_status%message)
        return
      end if
      valuation%policy_value = valuation%policy_value + one%policy_value
      valuation%risk_charge = valuation%risk_charge + one%risk_charge
    end do
    valuation%policy_value = valuation%policy_value / real(size(fund_scenarios,1),dp)
    valuation%risk_charge = valuation%risk_charge / real(size(fund_scenarios,1),dp)
  end subroutine valuate_one_policy_3d

  subroutine valuate_portfolio(portfolio, mortality_table, fund_scenarios, dt, discount, valuation, options, status)
    type(portfolio_type), intent(in) :: portfolio
    type(mortality_table_type), intent(in) :: mortality_table
    real(dp), intent(in) :: fund_scenarios(:,:,:), dt, discount(:)
    type(portfolio_valuation_type), intent(out) :: valuation
    type(valuation_options_type), intent(in), optional :: options
    type(status_type), intent(inout), optional :: status
    type(policy_valuation_type) :: one
    type(status_type) :: local_status
    integer :: i, n
    if (present(status)) call status%clear()
    n = portfolio%size()
    if (n < 1) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Portfolio is empty.')
      return
    end if
    allocate(valuation%policy_values(n), valuation%policy_risk_charges(n))
    do i = 1, n
      call valuate_one_policy_3d(portfolio%policies(i), mortality_table, fund_scenarios, dt, discount, one, options, local_status)
      if (.not. local_status%ok()) then
        if (present(status)) call status%set(local_status%code, local_status%message)
        return
      end if
      valuation%policy_values(i) = one%policy_value
      valuation%policy_risk_charges(i) = one%risk_charge
    end do
    valuation%portfolio_value = sum(valuation%policy_values)
    valuation%portfolio_risk_charge = sum(valuation%policy_risk_charges)
  end subroutine valuate_portfolio

  subroutine age_one_policy(policy, mortality_table, historical_fund_scenario, scenario_dates, dt, target_date, discount, &
                            aged_policy, options, status)
    type(policy_type), intent(in) :: policy
    type(mortality_table_type), intent(in) :: mortality_table
    real(dp), intent(in) :: historical_fund_scenario(:,:), dt, discount(:)
    type(date_type), intent(in) :: scenario_dates(:), target_date
    type(policy_type), intent(out) :: aged_policy
    type(valuation_options_type), intent(in), optional :: options
    type(status_type), intent(inout), optional :: status
    type(valuation_options_type) :: opt
    type(mortality_factors_type) :: mortality
    type(projection_result_type) :: projection
    type(status_type) :: local_status
    type(date_type) :: rolling_end
    integer :: start_index, end_index, nstep
    if (present(status)) call status%clear()
    opt = valuation_options_type()
    if (present(options)) opt = options
    aged_policy = policy
    if (target_date <= policy%current_date) return
    if (size(scenario_dates) < 1 .or. size(historical_fund_scenario,1) /= size(scenario_dates)) then
      if (present(status)) call status%set(vamc_dimension_error, 'Historical scenario dates and rows must align.')
      return
    end if
    if (target_date > scenario_dates(size(scenario_dates))) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Target date is beyond the last historical scenario date.')
      return
    end if
    if (policy%current_date < scenario_dates(1)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Current date is before the first historical scenario date.')
      return
    end if
    rolling_end = add_months(target_date,-1)
    if (opt%source_compatible_aging_index) then
      start_index = max(1, months_between(scenario_dates(1), policy%current_date))
      end_index = max(start_index, months_between(scenario_dates(1), rolling_end))
    else
      start_index = months_between(scenario_dates(1), policy%current_date) + 1
      end_index = months_between(scenario_dates(1), rolling_end) + 1
    end if
    start_index = min(max(start_index,1),size(scenario_dates))
    end_index = min(max(end_index,start_index),size(scenario_dates))
    nstep = end_index-start_index+1
    if (size(discount) < nstep) then
      if (present(status)) call status%set(vamc_dimension_error, 'Discount horizon is too short for aging.')
      return
    end if
    call calc_mort_factors(policy, mortality_table, dt, mortality, local_status)
    if (.not. local_status%ok() .or. size(mortality%pq) < nstep) then
      if (present(status)) call status%set(vamc_dimension_error, 'Mortality horizon is too short for aging.')
      return
    end if
    call project_policy(policy, historical_fund_scenario(start_index:end_index,:), dt, mortality%pq, mortality%p, discount, &
                        projection, opt, local_status)
    if (.not. local_status%ok()) then
      if (present(status)) call status%set(local_status%code, local_status%message)
      return
    end if
    aged_policy = projection%out_policy
  end subroutine age_one_policy

  subroutine age_portfolio(portfolio, mortality_table, historical_fund_scenario, scenario_dates, dt, target_date, discount, &
                           aged_portfolio, options, status)
    type(portfolio_type), intent(in) :: portfolio
    type(mortality_table_type), intent(in) :: mortality_table
    real(dp), intent(in) :: historical_fund_scenario(:,:), dt, discount(:)
    type(date_type), intent(in) :: scenario_dates(:), target_date
    type(portfolio_type), intent(out) :: aged_portfolio
    type(valuation_options_type), intent(in), optional :: options
    type(status_type), intent(inout), optional :: status
    type(status_type) :: local_status
    integer :: i, n
    if (present(status)) call status%clear()
    n = portfolio%size()
    allocate(aged_portfolio%policies(n))
    do i = 1, n
      call age_one_policy(portfolio%policies(i), mortality_table, historical_fund_scenario, scenario_dates, dt, target_date, &
                          discount, aged_portfolio%policies(i), options, local_status)
      if (.not. local_status%ok()) then
        if (present(status)) call status%set(local_status%code, local_status%message)
        return
      end if
    end do
  end subroutine age_portfolio
end module vamc_valuation
