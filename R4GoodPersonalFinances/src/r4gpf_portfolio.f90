! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_portfolio
  use r4gpf_kinds, only: dp, i8
  use r4gpf_status, only: r4gpf_success, r4gpf_invalid_argument, r4gpf_dimension_error, r4gpf_not_converged, &
    r4gpf_numerical_error
  use r4gpf_linalg, only: covariance_from_sd_corr, quadratic_form, project_simplex, cholesky_lower, vector_norm2
  use r4gpf_random, only: rng_state, seed_rng, random_normal
  use r4gpf_finance, only: tax_assumptions, tax_result, effective_tax_rates
  implicit none
  private
  public :: portfolio_spec, portfolio_result, portfolio_parameters
  public :: create_default_portfolio, calculate_portfolio_parameters, expected_utility
  public :: optimize_portfolio, networth_fractions, generate_random_returns
  public :: calculate_joint_networth_moments

  type :: portfolio_spec
    integer :: n_assets = 0
    character(len=64), allocatable :: names(:)
    real(dp), allocatable :: expected_return(:)
    real(dp), allocatable :: standard_deviation(:)
    real(dp), allocatable :: correlations(:, :)
    real(dp), allocatable :: taxable_accounts(:)
    real(dp), allocatable :: taxadvantaged_accounts(:)
    real(dp), allocatable :: human_capital_weights(:)
    real(dp), allocatable :: liabilities_weights(:)
    real(dp), allocatable :: effective_tax_rate(:)
    type(tax_assumptions) :: pretax
  end type portfolio_spec

  type :: portfolio_parameters
    real(dp) :: value = 0.0_dp
    real(dp), allocatable :: weights(:)
    real(dp) :: expected_return = 0.0_dp
    real(dp) :: standard_deviation = 0.0_dp
    integer :: status = r4gpf_success
  end type portfolio_parameters

  type :: portfolio_result
    real(dp), allocatable :: taxable(:)
    real(dp), allocatable :: taxadvantaged(:)
    real(dp), allocatable :: total(:)
    real(dp) :: expected_return = 0.0_dp
    real(dp) :: variance = 0.0_dp
    real(dp) :: standard_deviation = 0.0_dp
    real(dp) :: objective = -huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = r4gpf_success
  end type portfolio_result
contains

  subroutine create_default_portfolio(portfolio)
    type(portfolio_spec), intent(out) :: portfolio
    type(tax_result) :: taxes
    integer :: n

    n = 2
    portfolio%n_assets = n
    allocate(portfolio%names(n), portfolio%expected_return(n), portfolio%standard_deviation(n), &
      portfolio%correlations(n, n), portfolio%taxable_accounts(n), portfolio%taxadvantaged_accounts(n), &
      portfolio%human_capital_weights(n), portfolio%liabilities_weights(n), portfolio%effective_tax_rate(n))
    portfolio%names = [character(len=64) :: "GlobalStocksIndexFund", "InflationProtectedBonds"]
    portfolio%expected_return = [0.0461_dp, 0.0200_dp]
    portfolio%standard_deviation = [0.15_dp, 0.0_dp]
    portfolio%correlations = 0.0_dp
    portfolio%correlations(1, 1) = 1.0_dp
    portfolio%correlations(2, 2) = 1.0_dp
    portfolio%taxable_accounts = 0.0_dp
    portfolio%taxadvantaged_accounts = 0.0_dp
    portfolio%human_capital_weights = 0.5_dp
    portfolio%liabilities_weights = 0.5_dp
    allocate(portfolio%pretax%turnover(n), portfolio%pretax%income_qualified(n), &
      portfolio%pretax%capital_gains_long_term(n), portfolio%pretax%income(n), &
      portfolio%pretax%capital_gains(n), portfolio%pretax%cost_basis(n))
    portfolio%pretax%turnover = 0.04_dp
    portfolio%pretax%income_qualified = 0.0_dp
    portfolio%pretax%capital_gains_long_term = 1.0_dp
    portfolio%pretax%income = 0.0_dp
    portfolio%pretax%capital_gains = portfolio%expected_return
    portfolio%pretax%cost_basis = 1.0_dp / (1.0_dp + portfolio%expected_return)**10
    call effective_tax_rates(portfolio%expected_return, portfolio%standard_deviation, portfolio%pretax, &
      0.20_dp, 0.40_dp, taxes)
    portfolio%effective_tax_rate = taxes%effective_tax_rate
  end subroutine create_default_portfolio

  subroutine calculate_portfolio_parameters(portfolio, parameters)
    type(portfolio_spec), intent(in) :: portfolio
    type(portfolio_parameters), intent(out) :: parameters
    real(dp), allocatable :: covariance(:, :)
    real(dp) :: total_value
    integer :: st

    if (.not. validate_portfolio(portfolio)) then
      parameters%status = r4gpf_dimension_error
      return
    end if
    total_value = sum(portfolio%taxable_accounts) + sum(portfolio%taxadvantaged_accounts)
    if (total_value <= 0.0_dp) then
      parameters%status = r4gpf_invalid_argument
      return
    end if
    parameters%value = total_value
    allocate(parameters%weights(portfolio%n_assets))
    parameters%weights = (portfolio%taxable_accounts + portfolio%taxadvantaged_accounts) / total_value
    parameters%expected_return = dot_product(parameters%weights, portfolio%expected_return)
    call covariance_from_sd_corr(portfolio%standard_deviation, portfolio%correlations, covariance, st)
    if (st /= r4gpf_success) then
      parameters%status = st
      return
    end if
    parameters%standard_deviation = sqrt(max(0.0_dp, quadratic_form(parameters%weights, covariance)))
    parameters%status = r4gpf_success
  end subroutine calculate_portfolio_parameters

  elemental real(dp) function expected_utility(expected_return_value, variance, risk_tolerance) result(value)
    real(dp), intent(in) :: expected_return_value, variance, risk_tolerance
    real(dp) :: base, u, u2
    base = max(1.0e-8_dp, 1.0_dp + expected_return_value)
    if (risk_tolerance <= 0.0_dp) then
      value = -huge(1.0_dp)
      return
    end if
    if (abs(risk_tolerance - 1.0_dp) <= 10.0_dp * epsilon(1.0_dp)) then
      u = log(base)
    else
      u = risk_tolerance / (risk_tolerance - 1.0_dp) * base**((risk_tolerance - 1.0_dp) / risk_tolerance)
    end if
    u2 = -1.0_dp / (risk_tolerance * base**((1.0_dp + risk_tolerance) / risk_tolerance))
    value = u + 0.5_dp * u2 * variance
  end function expected_utility

  subroutine networth_fractions(financial_wealth, human_capital, liabilities, nondiscretionary_consumption, &
      discretionary_consumption, income, life_insurance_premium, financial_fraction, human_fraction, &
      liabilities_fraction, status)
    real(dp), intent(in) :: financial_wealth, human_capital, liabilities, nondiscretionary_consumption
    real(dp), intent(in) :: discretionary_consumption, income, life_insurance_premium
    real(dp), intent(out) :: financial_fraction, human_fraction, liabilities_fraction
    integer, intent(out) :: status
    real(dp) :: fw_prime, hc_prime, liab_prime, nw_prime

    fw_prime = financial_wealth + income - discretionary_consumption - nondiscretionary_consumption - &
      life_insurance_premium
    hc_prime = human_capital - income
    liab_prime = liabilities - nondiscretionary_consumption - life_insurance_premium
    nw_prime = fw_prime + hc_prime - liab_prime
    if (abs(nw_prime) <= 100.0_dp * epsilon(1.0_dp)) then
      financial_fraction = 0.0_dp
      human_fraction = 0.0_dp
      liabilities_fraction = 0.0_dp
      status = r4gpf_numerical_error
      return
    end if
    financial_fraction = fw_prime / nw_prime
    human_fraction = hc_prime / nw_prime
    liabilities_fraction = liab_prime / nw_prime
    status = r4gpf_success
  end subroutine networth_fractions

  subroutine optimize_portfolio(risk_tolerance, expected_returns, standard_deviations, correlations, result, &
      effective_tax_rates_input, fraction_taxable, financial_wealth, human_capital, liabilities, &
      nondiscretionary_consumption, discretionary_consumption, income, life_insurance_premium, &
      human_capital_weights, liabilities_weights, initial_allocation, max_iterations, tolerance)
    real(dp), intent(in) :: risk_tolerance, expected_returns(:), standard_deviations(:), correlations(:, :)
    type(portfolio_result), intent(out) :: result
    real(dp), intent(in), optional :: effective_tax_rates_input(:), fraction_taxable
    real(dp), intent(in), optional :: financial_wealth, human_capital, liabilities, nondiscretionary_consumption
    real(dp), intent(in), optional :: discretionary_consumption, income, life_insurance_premium
    real(dp), intent(in), optional :: human_capital_weights(:), liabilities_weights(:), initial_allocation(:)
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: covariance(:, :), w(:), gradient(:), candidate(:), direction(:), previous_w(:)
    real(dp), allocatable :: tax_multiplier(:), y(:), grad_mu(:), grad_var(:), cvec(:)
    real(dp) :: target_tax, target_adv, mu, var, obj, candidate_obj, step, tol
    real(dp) :: fw_frac, hc_frac, liab_frac, u1, u2, u3, base
    integer :: n, nw, iter, maxit, st, line_iter
    logical :: has_tax, has_networth

    n = size(expected_returns)
    if (n < 1 .or. size(standard_deviations) /= n .or. size(correlations, 1) /= n .or. size(correlations, 2) /= n) then
      result%status = r4gpf_dimension_error
      return
    end if
    if (risk_tolerance <= 0.0_dp) then
      result%status = r4gpf_invalid_argument
      return
    end if
    has_tax = present(effective_tax_rates_input)
    has_networth = present(financial_wealth)
    if (has_tax .neqv. present(fraction_taxable)) then
      result%status = r4gpf_invalid_argument
      return
    end if
    if (has_networth) then
      if (.not. has_tax .or. .not. present(human_capital) .or. .not. present(liabilities) .or. &
          .not. present(nondiscretionary_consumption) .or. .not. present(discretionary_consumption) .or. &
          .not. present(income) .or. .not. present(life_insurance_premium) .or. &
          .not. present(human_capital_weights) .or. .not. present(liabilities_weights)) then
        result%status = r4gpf_invalid_argument
        return
      end if
      if (size(human_capital_weights) /= n .or. size(liabilities_weights) /= n) then
        result%status = r4gpf_dimension_error
        return
      end if
    end if
    if (has_tax) then
      if (size(effective_tax_rates_input) /= n) then
        result%status = r4gpf_dimension_error
        return
      end if
      nw = 2 * n
      target_tax = min(1.0_dp, max(0.0_dp, fraction_taxable))
      target_adv = 1.0_dp - target_tax
      allocate(tax_multiplier(n))
      tax_multiplier = 1.0_dp - effective_tax_rates_input
    else
      nw = n
      target_tax = 1.0_dp
      target_adv = 0.0_dp
      allocate(tax_multiplier(n))
      tax_multiplier = 1.0_dp
    end if
    call covariance_from_sd_corr(standard_deviations, correlations, covariance, st)
    if (st /= r4gpf_success) then
      result%status = st
      return
    end if
    allocate(w(nw), gradient(nw), candidate(nw), direction(nw), previous_w(nw), y(n), grad_mu(nw), grad_var(nw), cvec(n))
    if (present(initial_allocation)) then
      if (size(initial_allocation) == nw) then
        w = initial_allocation
      else
        w = 1.0_dp / real(nw, dp)
      end if
    else
      w = 1.0_dp / real(nw, dp)
    end if
    call project_blocks(w)
    fw_frac = 1.0_dp
    hc_frac = 0.0_dp
    liab_frac = 0.0_dp
    cvec = 0.0_dp
    if (has_networth) then
      call networth_fractions(financial_wealth, human_capital, liabilities, nondiscretionary_consumption, &
        discretionary_consumption, income, life_insurance_premium, fw_frac, hc_frac, liab_frac, st)
      if (st /= r4gpf_success) then
        result%status = st
        return
      end if
      cvec = hc_frac * matmul(covariance, human_capital_weights) - &
        liab_frac * matmul(covariance, liabilities_weights)
    end if
    tol = 1.0e-11_dp
    if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
    maxit = 20000
    if (present(max_iterations)) maxit = max(10, max_iterations)
    call evaluate(w, mu, var, obj, gradient)
    step = 1.0_dp
    do iter = 1, maxit
      previous_w = w
      candidate = w + step * gradient
      call project_blocks(candidate)
      direction = candidate - w
      if (vector_norm2(direction) <= tol * (1.0_dp + vector_norm2(w))) exit
      candidate_obj = objective_only(candidate)
      line_iter = 0
      do while (candidate_obj < obj + 1.0e-4_dp * dot_product(gradient, direction) .and. line_iter < 50)
        step = step * 0.5_dp
        candidate = w + step * gradient
        call project_blocks(candidate)
        direction = candidate - w
        candidate_obj = objective_only(candidate)
        line_iter = line_iter + 1
      end do
      if (candidate_obj <= obj + 10.0_dp * epsilon(1.0_dp)) then
        exit
      end if
      w = candidate
      call evaluate(w, mu, var, obj, gradient)
      if (vector_norm2(w - previous_w) <= tol * (1.0_dp + vector_norm2(w))) exit
      step = min(100.0_dp, step * 1.25_dp)
    end do
    allocate(result%taxable(n), result%taxadvantaged(n), result%total(n))
    if (has_tax) then
      result%taxable = w(1:n)
      result%taxadvantaged = w(n + 1:2 * n)
    else
      result%taxable = w
      result%taxadvantaged = 0.0_dp
    end if
    result%total = result%taxable + result%taxadvantaged
    result%expected_return = mu
    result%variance = max(0.0_dp, var)
    result%standard_deviation = sqrt(result%variance)
    result%objective = obj
    result%iterations = min(iter, maxit)
    if (iter <= maxit) then
      result%status = r4gpf_success
    else
      result%status = r4gpf_not_converged
    end if

  contains
    subroutine project_blocks(x)
      real(dp), intent(inout) :: x(:)
      real(dp), allocatable :: p(:)
      allocate(p(n))
      if (has_tax) then
        call project_simplex(x(1:n), target_tax, p)
        x(1:n) = p
        call project_simplex(x(n + 1:2 * n), target_adv, p)
        x(n + 1:2 * n) = p
      else
        call project_simplex(x, 1.0_dp, p)
        x = p
      end if
    end subroutine project_blocks

    subroutine moments(x, mean_value, variance_value, mean_gradient, variance_gradient)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: mean_value, variance_value
      real(dp), intent(out) :: mean_gradient(:), variance_gradient(:)
      real(dp) :: constant_mean, constant_var
      real(dp), allocatable :: gv_y(:)
      allocate(gv_y(n))
      if (has_tax) then
        y = tax_multiplier * x(1:n) + x(n + 1:2 * n)
      else
        y = x
      end if
      constant_mean = 0.0_dp
      constant_var = 0.0_dp
      if (has_networth) then
        constant_mean = hc_frac * dot_product(human_capital_weights, expected_returns) - &
          liab_frac * dot_product(liabilities_weights, expected_returns)
        constant_var = hc_frac**2 * quadratic_form(human_capital_weights, covariance) + &
          liab_frac**2 * quadratic_form(liabilities_weights, covariance) - &
          2.0_dp * hc_frac * liab_frac * dot_product(human_capital_weights, matmul(covariance, liabilities_weights))
        mean_value = fw_frac * dot_product(y, expected_returns) + constant_mean
        variance_value = fw_frac**2 * quadratic_form(y, covariance) + &
          2.0_dp * fw_frac * dot_product(cvec, y) + constant_var
        gv_y = 2.0_dp * fw_frac**2 * matmul(covariance, y) + 2.0_dp * fw_frac * cvec
        if (has_tax) then
          mean_gradient(1:n) = fw_frac * tax_multiplier * expected_returns
          mean_gradient(n + 1:2 * n) = fw_frac * expected_returns
          variance_gradient(1:n) = tax_multiplier * gv_y
          variance_gradient(n + 1:2 * n) = gv_y
        else
          mean_gradient = fw_frac * expected_returns
          variance_gradient = gv_y
        end if
      else
        mean_value = dot_product(y, expected_returns)
        variance_value = quadratic_form(y, covariance)
        gv_y = 2.0_dp * matmul(covariance, y)
        if (has_tax) then
          mean_gradient(1:n) = tax_multiplier * expected_returns
          mean_gradient(n + 1:2 * n) = expected_returns
          variance_gradient(1:n) = tax_multiplier * gv_y
          variance_gradient(n + 1:2 * n) = gv_y
        else
          mean_gradient = expected_returns
          variance_gradient = gv_y
        end if
      end if
    end subroutine moments

    subroutine evaluate(x, mean_value, variance_value, objective_value, objective_gradient)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: mean_value, variance_value, objective_value
      real(dp), intent(out) :: objective_gradient(:)
      call moments(x, mean_value, variance_value, grad_mu, grad_var)
      base = max(1.0e-8_dp, 1.0_dp + mean_value)
      u1 = base**(-1.0_dp / risk_tolerance)
      u2 = -1.0_dp / (risk_tolerance * base**((1.0_dp + risk_tolerance) / risk_tolerance))
      u3 = (1.0_dp + risk_tolerance) / risk_tolerance**2 * &
        base**(-(1.0_dp + 2.0_dp * risk_tolerance) / risk_tolerance)
      objective_value = expected_utility(mean_value, variance_value, risk_tolerance)
      objective_gradient = (u1 + 0.5_dp * u3 * variance_value) * grad_mu + 0.5_dp * u2 * grad_var
    end subroutine evaluate

    function objective_only(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, mm, vv
      call moments(x, mm, vv, grad_mu, grad_var)
      value = expected_utility(mm, vv, risk_tolerance)
    end function objective_only
  end subroutine optimize_portfolio

  subroutine generate_random_returns(portfolio, n_periods, returns, status, seed)
    type(portfolio_spec), intent(in) :: portfolio
    integer, intent(in) :: n_periods
    real(dp), allocatable, intent(out) :: returns(:, :)
    integer, intent(out) :: status
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: covariance(:, :), adjusted_sd(:), l(:, :), z(:)
    type(rng_state) :: rng
    integer :: i, j, st
    logical, allocatable :: zero_sd(:)

    if (.not. validate_portfolio(portfolio) .or. n_periods < 1) then
      allocate(returns(0, 0))
      status = r4gpf_invalid_argument
      return
    end if
    allocate(adjusted_sd(portfolio%n_assets), zero_sd(portfolio%n_assets), z(portfolio%n_assets))
    adjusted_sd = portfolio%standard_deviation
    zero_sd = adjusted_sd <= 0.0_dp
    where (zero_sd) adjusted_sd = 1.0e-7_dp
    call covariance_from_sd_corr(adjusted_sd, portfolio%correlations, covariance, st)
    if (st /= r4gpf_success) then
      allocate(returns(0, 0))
      status = st
      return
    end if
    call cholesky_lower(covariance, l, st, 1.0e-14_dp)
    if (st /= r4gpf_success) then
      allocate(returns(0, 0))
      status = st
      return
    end if
    if (present(seed)) then
      call seed_rng(rng, seed)
    else
      call seed_rng(rng, 123456789_i8)
    end if
    allocate(returns(n_periods, portfolio%n_assets))
    do i = 1, n_periods
      do j = 1, portfolio%n_assets
        z(j) = random_normal(rng)
      end do
      returns(i, :) = portfolio%expected_return + matmul(l, z)
      where (zero_sd) returns(i, :) = portfolio%expected_return
    end do
    status = r4gpf_success
  end subroutine generate_random_returns


  subroutine calculate_joint_networth_moments(allocations_taxable, allocations_taxadvantaged, expected_returns, &
      standard_deviations, correlations, effective_tax_rates_input, human_capital_weights, liabilities_weights, &
      financial_wealth, human_capital, liabilities, nondiscretionary_consumption, discretionary_consumption, &
      income, life_insurance_premium, mean_value, variance_value, status)
    real(dp), intent(in) :: allocations_taxable(:), allocations_taxadvantaged(:), expected_returns(:)
    real(dp), intent(in) :: standard_deviations(:), correlations(:, :), effective_tax_rates_input(:)
    real(dp), intent(in) :: human_capital_weights(:), liabilities_weights(:)
    real(dp), intent(in) :: financial_wealth, human_capital, liabilities, nondiscretionary_consumption
    real(dp), intent(in) :: discretionary_consumption, income, life_insurance_premium
    real(dp), intent(out) :: mean_value, variance_value
    integer, intent(out) :: status
    real(dp), allocatable :: covariance(:, :), y(:), cvec(:)
    real(dp) :: fw_frac, hc_frac, liab_frac
    integer :: n, st

    n = size(expected_returns)
    if (size(allocations_taxable) /= n .or. size(allocations_taxadvantaged) /= n .or. &
        size(standard_deviations) /= n .or. size(effective_tax_rates_input) /= n .or. &
        size(human_capital_weights) /= n .or. size(liabilities_weights) /= n .or. &
        size(correlations, 1) /= n .or. size(correlations, 2) /= n) then
      mean_value = 0.0_dp
      variance_value = 0.0_dp
      status = r4gpf_dimension_error
      return
    end if
    call covariance_from_sd_corr(standard_deviations, correlations, covariance, st)
    if (st /= r4gpf_success) then
      mean_value = 0.0_dp
      variance_value = 0.0_dp
      status = st
      return
    end if
    call networth_fractions(financial_wealth, human_capital, liabilities, nondiscretionary_consumption, &
      discretionary_consumption, income, life_insurance_premium, fw_frac, hc_frac, liab_frac, st)
    if (st /= r4gpf_success) then
      mean_value = 0.0_dp
      variance_value = 0.0_dp
      status = st
      return
    end if
    allocate(y(n), cvec(n))
    y = (1.0_dp - effective_tax_rates_input) * allocations_taxable + allocations_taxadvantaged
    cvec = hc_frac * matmul(covariance, human_capital_weights) - &
      liab_frac * matmul(covariance, liabilities_weights)
    mean_value = fw_frac * dot_product(y, expected_returns) + &
      hc_frac * dot_product(human_capital_weights, expected_returns) - &
      liab_frac * dot_product(liabilities_weights, expected_returns)
    variance_value = fw_frac**2 * quadratic_form(y, covariance) + 2.0_dp * fw_frac * dot_product(cvec, y) + &
      hc_frac**2 * quadratic_form(human_capital_weights, covariance) + &
      liab_frac**2 * quadratic_form(liabilities_weights, covariance) - &
      2.0_dp * hc_frac * liab_frac * dot_product(human_capital_weights, matmul(covariance, liabilities_weights))
    status = r4gpf_success
  end subroutine calculate_joint_networth_moments

  logical function validate_portfolio(portfolio) result(valid)
    type(portfolio_spec), intent(in) :: portfolio
    integer :: n
    n = portfolio%n_assets
    valid = n > 0 .and. allocated(portfolio%expected_return) .and. allocated(portfolio%standard_deviation) .and. &
      allocated(portfolio%correlations) .and. allocated(portfolio%taxable_accounts) .and. &
      allocated(portfolio%taxadvantaged_accounts) .and. allocated(portfolio%human_capital_weights) .and. &
      allocated(portfolio%liabilities_weights) .and. allocated(portfolio%effective_tax_rate)
    if (.not. valid) return
    valid = size(portfolio%expected_return) == n .and. size(portfolio%standard_deviation) == n .and. &
      size(portfolio%correlations, 1) == n .and. size(portfolio%correlations, 2) == n .and. &
      size(portfolio%taxable_accounts) == n .and. size(portfolio%taxadvantaged_accounts) == n .and. &
      size(portfolio%human_capital_weights) == n .and. size(portfolio%liabilities_weights) == n .and. &
      size(portfolio%effective_tax_rate) == n
  end function validate_portfolio

end module r4gpf_portfolio
