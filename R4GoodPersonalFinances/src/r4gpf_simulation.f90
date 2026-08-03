! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_simulation
  use r4gpf_kinds, only: dp, i8
  use r4gpf_status, only: r4gpf_success, r4gpf_invalid_argument, r4gpf_dimension_error
  use r4gpf_finance, only: present_value_stream, certainty_equivalent_return, utility
  use r4gpf_mortality, only: gompertz_survival_probability
  use r4gpf_household, only: household, date_type, household_timeline, build_household_timeline
  use r4gpf_portfolio, only: portfolio_spec, portfolio_result, optimize_portfolio, &
    calculate_joint_networth_moments, generate_random_returns
  implicit none
  private
  public :: lifecycle_result, discretionary_spending, consumption_growth_rate, consumption_delta
  public :: simulate_lifecycle, simulate_lifecycle_samples

  type :: lifecycle_result
    integer :: n_periods = 0
    integer :: n_assets = 0
    integer :: status = r4gpf_success
    type(household_timeline) :: timeline
    real(dp), allocatable :: income(:)
    real(dp), allocatable :: nondiscretionary_spending(:)
    real(dp), allocatable :: human_capital(:)
    real(dp), allocatable :: liabilities(:)
    real(dp), allocatable :: financial_wealth(:)
    real(dp), allocatable :: net_worth(:)
    real(dp), allocatable :: discretionary_spending(:)
    real(dp), allocatable :: total_spending(:)
    real(dp), allocatable :: financial_wealth_end(:)
    real(dp), allocatable :: returns(:, :)
    real(dp), allocatable :: taxable_allocation(:, :)
    real(dp), allocatable :: taxadvantaged_allocation(:, :)
    real(dp), allocatable :: total_allocation(:, :)
    real(dp), allocatable :: time_value_discount(:)
    real(dp), allocatable :: spending_utility(:)
    real(dp), allocatable :: weighted_spending_utility(:)
    real(dp), allocatable :: savings(:)
    real(dp), allocatable :: saving_rate(:)
  end type lifecycle_result
contains

  elemental real(dp) function consumption_growth_rate(discount_rate, consumption_impatience_preference, &
      smooth_consumption_preference) result(value)
    real(dp), intent(in) :: discount_rate, consumption_impatience_preference, smooth_consumption_preference
    if (1.0_dp + consumption_impatience_preference <= 0.0_dp) then
      value = 0.0_dp
    else
      value = ((1.0_dp + discount_rate) / (1.0_dp + consumption_impatience_preference))** &
        smooth_consumption_preference - 1.0_dp
    end if
  end function consumption_growth_rate

  real(dp) function consumption_delta(survival_probabilities, smooth_consumption_preference, growth_rate, &
      discount_rate) result(value)
    real(dp), intent(in) :: survival_probabilities(:), smooth_consumption_preference, growth_rate, discount_rate
    real(dp) :: ratio
    integer :: i
    ratio = (1.0_dp + growth_rate) / (1.0_dp + discount_rate)
    if (.not. (ratio > -huge(1.0_dp) .and. ratio < huge(1.0_dp))) ratio = 1.0_dp
    value = 0.0_dp
    do i = 1, size(survival_probabilities)
      value = value + survival_probabilities(i)**smooth_consumption_preference * ratio**real(i - 1, dp)
    end do
  end function consumption_delta

  real(dp) function discretionary_spending(net_worth, consumption_impatience_preference, &
      smooth_consumption_preference, current_age, max_age, gompertz_mode, gompertz_dispersion, &
      discount_rate, allocations_taxable, allocations_taxadvantaged, expected_returns, standard_deviations, &
      effective_tax_rates, correlations, financial_wealth, human_capital, human_capital_weights, liabilities, &
      liabilities_weights, income, nondiscretionary_consumption, life_insurance_premium, risk_tolerance, status) &
      result(value)
    real(dp), intent(in) :: net_worth, consumption_impatience_preference, smooth_consumption_preference
    real(dp), intent(in) :: current_age, max_age, gompertz_mode, gompertz_dispersion
    real(dp), intent(in), optional :: discount_rate
    real(dp), intent(in), optional :: allocations_taxable(:), allocations_taxadvantaged(:), expected_returns(:)
    real(dp), intent(in), optional :: standard_deviations(:), effective_tax_rates(:), correlations(:, :)
    real(dp), intent(in), optional :: financial_wealth, human_capital, human_capital_weights(:), liabilities
    real(dp), intent(in), optional :: liabilities_weights(:), income, nondiscretionary_consumption
    real(dp), intent(in), optional :: life_insurance_premium, risk_tolerance
    integer, intent(out), optional :: status
    real(dp), allocatable :: survival(:)
    real(dp) :: rate, mean_value, variance_value, growth, delta
    integer :: i, n, st

    if (present(discount_rate)) then
      rate = discount_rate
      st = r4gpf_success
    else
      if (.not. present(allocations_taxable) .or. .not. present(allocations_taxadvantaged) .or. &
          .not. present(expected_returns) .or. .not. present(standard_deviations) .or. &
          .not. present(effective_tax_rates) .or. .not. present(correlations) .or. &
          .not. present(financial_wealth) .or. .not. present(human_capital) .or. &
          .not. present(human_capital_weights) .or. .not. present(liabilities) .or. &
          .not. present(liabilities_weights) .or. .not. present(income) .or. &
          .not. present(nondiscretionary_consumption) .or. .not. present(life_insurance_premium) .or. &
          .not. present(risk_tolerance)) then
        value = 0.0_dp
        if (present(status)) status = r4gpf_invalid_argument
        return
      end if
      call calculate_joint_networth_moments(allocations_taxable, allocations_taxadvantaged, expected_returns, &
        standard_deviations, correlations, effective_tax_rates, human_capital_weights, liabilities_weights, &
        financial_wealth, human_capital, liabilities, nondiscretionary_consumption, 0.0_dp, income, &
        life_insurance_premium, mean_value, variance_value, st)
      if (st /= r4gpf_success) then
        value = 0.0_dp
        if (present(status)) status = st
        return
      end if
      rate = certainty_equivalent_return(mean_value, variance_value, risk_tolerance)
    end if
    n = max(1, floor(max_age - current_age) + 1)
    allocate(survival(n))
    do i = 1, n
      survival(i) = gompertz_survival_probability(current_age, current_age + real(i - 1, dp), &
        gompertz_mode, gompertz_dispersion, max_age)
    end do
    growth = consumption_growth_rate(rate, consumption_impatience_preference, smooth_consumption_preference)
    delta = consumption_delta(survival, smooth_consumption_preference, growth, rate)
    if (abs(delta) <= tiny(1.0_dp)) then
      value = net_worth
    else
      value = net_worth / delta
    end if
    if (present(status)) status = st
  end function discretionary_spending

  subroutine simulate_lifecycle(home, portfolio, current_date, income, nondiscretionary_spending, result, &
      use_random_returns, seed, optimizer_max_iterations)
    type(household), intent(in) :: home
    type(portfolio_spec), intent(in) :: portfolio
    type(date_type), intent(in) :: current_date
    real(dp), intent(in) :: income(:), nondiscretionary_spending(:)
    type(lifecycle_result), intent(out) :: result
    logical, intent(in), optional :: use_random_returns
    integer(i8), intent(in), optional :: seed
    integer, intent(in), optional :: optimizer_max_iterations
    real(dp), allocatable :: pv_income(:), pv_spending(:), initial_taxable(:), initial_advantaged(:), initial(:)
    real(dp) :: financial_wealth_initial, hc_rate, liab_rate, fraction_taxable, spending, wealth_before_return
    real(dp) :: min_age, max_age
    type(portfolio_result) :: optimal
    integer :: i, n, n_assets, st, maxit
    logical :: randomize

    call build_household_timeline(home, current_date, result%timeline, st)
    if (st /= r4gpf_success) then
      result%status = st
      return
    end if
    n = result%timeline%n_periods
    n_assets = portfolio%n_assets
    if (size(income) /= n .or. size(nondiscretionary_spending) /= n .or. n_assets < 1) then
      result%status = r4gpf_dimension_error
      return
    end if
    result%n_periods = n
    result%n_assets = n_assets
    allocate(result%income(n), result%nondiscretionary_spending(n), result%human_capital(n), result%liabilities(n), &
      result%financial_wealth(n), result%net_worth(n), result%discretionary_spending(n), result%total_spending(n), &
      result%financial_wealth_end(n), result%returns(n, n_assets), result%taxable_allocation(n_assets, n), &
      result%taxadvantaged_allocation(n_assets, n), result%total_allocation(n_assets, n), &
      result%time_value_discount(n), result%spending_utility(n), result%weighted_spending_utility(n), &
      result%savings(n), result%saving_rate(n), initial_taxable(n_assets), initial_advantaged(n_assets), &
      initial(2 * n_assets))
    result%income = income
    result%nondiscretionary_spending = nondiscretionary_spending
    financial_wealth_initial = sum(portfolio%taxable_accounts) + sum(portfolio%taxadvantaged_accounts)
    if (financial_wealth_initial > 0.0_dp) then
      initial_taxable = portfolio%taxable_accounts / financial_wealth_initial
      initial_advantaged = portfolio%taxadvantaged_accounts / financial_wealth_initial
      fraction_taxable = sum(portfolio%taxable_accounts) / financial_wealth_initial
    else
      initial_taxable = 0.0_dp
      initial_advantaged = 0.0_dp
      fraction_taxable = 0.0_dp
    end if
    hc_rate = dot_product(portfolio%human_capital_weights, portfolio%expected_return)
    liab_rate = dot_product(portfolio%liabilities_weights, portfolio%expected_return)
    call present_value_stream(income, hc_rate, pv_income, st)
    if (st /= r4gpf_success) then
      result%status = st
      return
    end if
    call present_value_stream(nondiscretionary_spending, liab_rate, pv_spending, st)
    if (st /= r4gpf_success) then
      result%status = st
      return
    end if
    result%human_capital = pv_income
    result%liabilities = pv_spending
    randomize = .false.
    if (present(use_random_returns)) randomize = use_random_returns
    if (randomize) then
      if (present(seed)) then
        call generate_random_returns(portfolio, n, result%returns, st, seed)
      else
        call generate_random_returns(portfolio, n, result%returns, st)
      end if
      if (st /= r4gpf_success) then
        result%status = st
        return
      end if
    else
      do i = 1, n
        result%returns(i, :) = portfolio%expected_return
      end do
    end if
    result%financial_wealth = 0.0_dp
    result%financial_wealth(1) = financial_wealth_initial
    result%taxable_allocation = 0.0_dp
    result%taxadvantaged_allocation = 0.0_dp
    result%total_allocation = 0.0_dp
    initial(1:n_assets) = initial_taxable
    initial(n_assets + 1:2 * n_assets) = initial_advantaged
    min_age = result%timeline%joint_fit%current_age
    max_age = min_age + real(n - 1, dp)
    maxit = 5000
    if (present(optimizer_max_iterations)) maxit = max(100, optimizer_max_iterations)

    do i = 1, n
      result%net_worth(i) = result%financial_wealth(i) + result%human_capital(i) - result%liabilities(i)
      spending = discretionary_spending(result%net_worth(i), home%consumption_impatience_preference, &
        home%smooth_consumption_preference, min_age + real(i - 1, dp), max_age, &
        result%timeline%joint_fit%mode, result%timeline%joint_fit%dispersion, &
        allocations_taxable=initial_taxable, allocations_taxadvantaged=initial_advantaged, &
        expected_returns=portfolio%expected_return, standard_deviations=portfolio%standard_deviation, &
        effective_tax_rates=portfolio%effective_tax_rate, correlations=portfolio%correlations, &
        financial_wealth=result%financial_wealth(i), human_capital=result%human_capital(i), &
        human_capital_weights=portfolio%human_capital_weights, liabilities=result%liabilities(i), &
        liabilities_weights=portfolio%liabilities_weights, income=income(i), &
        nondiscretionary_consumption=nondiscretionary_spending(i), life_insurance_premium=0.0_dp, &
        risk_tolerance=home%risk_tolerance, status=st)
      if (st /= r4gpf_success) then
        result%status = st
        return
      end if
      if (.not. (spending > -huge(1.0_dp) .and. spending < huge(1.0_dp))) spending = result%net_worth(i)
      result%discretionary_spending(i) = spending
      result%total_spending(i) = spending + nondiscretionary_spending(i)
      wealth_before_return = result%financial_wealth(i) + income(i) - result%total_spending(i)
      call optimize_portfolio(home%risk_tolerance, portfolio%expected_return, portfolio%standard_deviation, &
        portfolio%correlations, optimal, effective_tax_rates_input=portfolio%effective_tax_rate, &
        fraction_taxable=fraction_taxable, financial_wealth=result%financial_wealth(i), &
        human_capital=result%human_capital(i), liabilities=result%liabilities(i), &
        nondiscretionary_consumption=nondiscretionary_spending(i), discretionary_consumption=spending, &
        income=income(i), life_insurance_premium=0.0_dp, human_capital_weights=portfolio%human_capital_weights, &
        liabilities_weights=portfolio%liabilities_weights, initial_allocation=initial, max_iterations=maxit)
      if (optimal%status /= r4gpf_success) then
        if (i == 1) then
          optimal%taxable = initial_taxable
          optimal%taxadvantaged = initial_advantaged
          optimal%total = initial_taxable + initial_advantaged
        else
          optimal%taxable = result%taxable_allocation(:, i - 1)
          optimal%taxadvantaged = result%taxadvantaged_allocation(:, i - 1)
          optimal%total = result%total_allocation(:, i - 1)
        end if
      end if
      result%taxable_allocation(:, i) = optimal%taxable
      result%taxadvantaged_allocation(:, i) = optimal%taxadvantaged
      result%total_allocation(:, i) = optimal%total
      initial(1:n_assets) = optimal%taxable
      initial(n_assets + 1:2 * n_assets) = optimal%taxadvantaged
      result%financial_wealth_end(i) = wealth_before_return * &
        sum(optimal%total * (1.0_dp + result%returns(i, :)))
      if (i < n) result%financial_wealth(i + 1) = result%financial_wealth_end(i)
      result%time_value_discount(i) = 1.0_dp / &
        (1.0_dp + home%consumption_impatience_preference)**real(i - 1, dp)
      result%spending_utility(i) = utility(spending, home%smooth_consumption_preference)
      result%weighted_spending_utility(i) = result%timeline%gompertz_survival(i) * &
        result%time_value_discount(i) * result%spending_utility(i)
      if (income(i) > 0.0_dp .and. result%total_spending(i) > 0.0_dp) then
        result%savings(i) = max(0.0_dp, income(i) - result%total_spending(i))
        result%saving_rate(i) = max(0.0_dp, result%savings(i) / income(i))
      else
        result%savings(i) = 0.0_dp
        result%saving_rate(i) = 0.0_dp
      end if
    end do
    result%status = r4gpf_success
  end subroutine simulate_lifecycle


  subroutine simulate_lifecycle_samples(home, portfolio, current_date, income, nondiscretionary_spending, &
      n_samples, results, seeds, optimizer_max_iterations)
    type(household), intent(in) :: home
    type(portfolio_spec), intent(in) :: portfolio
    type(date_type), intent(in) :: current_date
    real(dp), intent(in) :: income(:), nondiscretionary_spending(:)
    integer, intent(in) :: n_samples
    type(lifecycle_result), allocatable, intent(out) :: results(:)
    integer(i8), intent(in), optional :: seeds(:)
    integer, intent(in), optional :: optimizer_max_iterations
    integer :: i
    integer(i8) :: sample_seed

    if (n_samples < 1) then
      allocate(results(0))
      return
    end if
    if (present(seeds)) then
      if (size(seeds) < n_samples) then
        allocate(results(0))
        return
      end if
    end if
    allocate(results(n_samples))
    do i = 1, n_samples
      if (present(seeds)) then
        sample_seed = seeds(i)
      else
        sample_seed = 104729_i8 * int(i, i8) + 12345_i8
      end if
      if (present(optimizer_max_iterations)) then
        call simulate_lifecycle(home, portfolio, current_date, income, nondiscretionary_spending, results(i), &
          use_random_returns=.true., seed=sample_seed, optimizer_max_iterations=optimizer_max_iterations)
      else
        call simulate_lifecycle(home, portfolio, current_date, income, nondiscretionary_spending, results(i), &
          use_random_returns=.true., seed=sample_seed)
      end if
    end do
  end subroutine simulate_lifecycle_samples

end module r4gpf_simulation
