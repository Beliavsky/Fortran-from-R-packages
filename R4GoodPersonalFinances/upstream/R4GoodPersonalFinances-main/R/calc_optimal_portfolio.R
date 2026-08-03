calc_optimal_portfolio <- function(
  risk_tolerance,
  expected_returns,
  standard_deviations,
  correlations,
  asset_names                  = NULL,
  effective_tax_rates          = NULL,
  in_taxable_accounts          = NULL,
  financial_wealth             = NULL,
  human_capital                = NULL,
  human_capital_weights        = NULL,
  liabilities                  = NULL,
  liabilities_weights          = NULL,
  nondiscretionary_consumption = NULL,
  discretionary_consumption    = NULL,
  income                       = NULL,
  life_insurance_premium       = NULL,
  initial_allocation           = NULL,
  opts = list(
    algorithm = "NLOPT_LD_SLSQP",
    maxeval   = 100
  )
) {

  taxable <- taxadvantaged <- asset_class <- NULL

  if (!is.null(in_taxable_accounts) && is.nan(in_taxable_accounts)) {
    in_taxable_accounts <- 0
  }

  covariance_matrix <- calc_covariance_matrix(
    standard_deviations = standard_deviations,
    correlations        = correlations
  )

  if (!is.null(effective_tax_rates)) {
    tax_matrix <- diag(1 - effective_tax_rates)
  }

  final_results <- new.env()
    
  objective_function <- function(params) {

    if (is.null(effective_tax_rates)) {

      expected_return <- calc_mvo_portfolio_expected_return(
        params           = params,
        expected_returns = expected_returns
      )
      variance <- calc_mvo_portfolio_variance(
        params            = params,
        covariance_matrix = covariance_matrix
      )
    }

    if (!is.null(effective_tax_rates) && is.null(financial_wealth)) {

      expected_return <- calc_joint_portfolio_expected_return(
        params           = params,
        expected_returns = expected_returns,
        tax_matrix       = tax_matrix
      )
      variance <- calc_joint_portfolio_variance(
        params            = params,
        covariance_matrix = covariance_matrix,
        tax_matrix        = tax_matrix
      )
    }
    
    if (!is.null(effective_tax_rates) && !is.null(financial_wealth)) {
          
      expected_return <- calc_joint_networth_portfolio_expected_return(
        params                       = params,
        expected_returns             = expected_returns,
        tax_matrix                   = tax_matrix,
        human_capital_weights        = human_capital_weights,
        liabilities_weights          = liabilities_weights,
        financial_wealth             = financial_wealth,
        human_capital                = human_capital,
        liabilities                  = liabilities,
        nondiscretionary_consumption = nondiscretionary_consumption,
        discretionary_consumption    = discretionary_consumption,
        income                       = income,
        life_insurance_premium       = life_insurance_premium
      )

      variance <- calc_joint_networth_portfolio_variance(
        params                       = params,
        covariance_matrix            = covariance_matrix,
        tax_matrix                   = tax_matrix,
        human_capital_weights        = human_capital_weights,
        liabilities_weights          = liabilities_weights,
        financial_wealth             = financial_wealth,
        human_capital                = human_capital,
        liabilities                  = liabilities,
        nondiscretionary_consumption = nondiscretionary_consumption,
        discretionary_consumption    = discretionary_consumption,
        income                       = income,
        life_insurance_premium       = life_insurance_premium
      )
    }

    final_results$expected_return <- expected_return
    final_results$variance        <- variance

    expected_utility <- calc_expected_utility(
      expected_return = expected_return,
      variance        = variance,
      risk_tolerance  = risk_tolerance
    )

    -expected_utility
  }
  
  assets_number <- length(expected_returns)
  
  if (is.null(effective_tax_rates)) {

    total_assets <- assets_number
    
    # Equality constraint: sum(allocations) = 1
    equality_constraint <- function(params) {
      return(sum(params) - 1)
    }

    equality_constraint_jacobian <- function(params) {
      n   <- length(params)
      jac <- matrix(1, nrow = 1, ncol = n)
      jac
    }

  } else {

    total_assets <- assets_number * 2 

    # Equality constraints:
    # 1. sum(allocations_taxable) = in_taxable_accounts
    # 2. sum(allocations_taxadvantaged) = 1 - in_taxable_accounts
    equality_constraint <- function(params) {

      allocations_taxable <- get_allocations_taxable(params)
      allocations_taxadvantaged <- get_allocations_taxadvantaged(params)

      equality_constraint <- c(
        sum(allocations_taxable) - in_taxable_accounts, 
        sum(allocations_taxadvantaged) - (1 - in_taxable_accounts)
      )

      equality_constraint
    }
    
    equality_constraint_jacobian <- function(params) {

      n   <- length(params)
      jac <- matrix(0, nrow = 2, ncol = n)
      jac[1, get_taxable_indices(params)] <- 1
      jac[2, get_taxadvantaged_indices(params)] <- 1
      jac
    }
  }

  if (is.null(initial_allocation)) {
    initial_allocation <- rep(1 / total_assets, total_assets)
  }

  # Set lower bounds for allocations (non-negativity)
  lower_bounds <- rep(0, total_assets)

  optimization_result <- nloptr::nloptr(
    x0          = initial_allocation,
    eval_f      = objective_function,
    eval_grad_f = function(initial_allocation) {
      nloptr::nl.grad(initial_allocation, objective_function)
    },
    opts          = opts,
    eval_g_eq     = equality_constraint,
    eval_jac_g_eq = equality_constraint_jacobian,
    lb            = lower_bounds
  )

  optimal_allocations <- optimization_result$solution

  if (is.null(effective_tax_rates)) {
    
    allocations <- 
      dplyr::tibble(
        total = optimal_allocations
      )

  } else {
    
    optimal_taxable_allocations <- 
      get_allocations_taxable(optimal_allocations)
    
    optimal_taxadvantaged_allocations <- 
      get_allocations_taxadvantaged(optimal_allocations)
    
    allocations <-
      dplyr::tibble(
        taxable       = optimal_taxable_allocations,
        taxadvantaged = optimal_taxadvantaged_allocations
      ) |>
        dplyr::mutate(
          total = taxable + taxadvantaged
        )
  }

  if (!is.null(asset_names)) {

    allocations <-
      allocations |> 
      dplyr::mutate(asset_class = asset_names) |> 
      dplyr::select(asset_class, dplyr::everything())
  }

  list(
    allocations        = allocations,
    expected_return    = as.numeric(final_results$expected_return),
    variance           = as.numeric(final_results$variance),
    standard_deviation = sqrt(as.numeric(final_results$variance))
  )
}

calc_expected_utility <- function(
  expected_return, 
  variance,
  risk_tolerance
) {

  if (risk_tolerance == 1) {
    crra_utility <- log(1 + expected_return)
  } else {

    # orginal formula
    # crra_utility <- 
    #   (risk_tolerance / (risk_tolerance - 1)) * 
    #   (1 + expected_return) ^ ((risk_tolerance - 1) / risk_tolerance)

    # stable gradient version
    crra_utility <- 
      (risk_tolerance / (risk_tolerance - 1)) * 
      (pmax(1 + expected_return, 1e-8)) ^ (
        (risk_tolerance - 1) / risk_tolerance
      )

  }
  
  # orginal formula
  # crra_utility_second_derivative <- 
  #   -1 / (
  #     risk_tolerance * 
  #       (1 + expected_return) ^ ((1 + risk_tolerance) / risk_tolerance)
  #   )
  
  # stable gradient version
  crra_utility_second_derivative <- 
    -1 / (
      risk_tolerance * 
        (pmax(1 + expected_return, 1e-8)) ^ (
          (1 + risk_tolerance) / risk_tolerance
        )
    )
  
  expected_utility <- 
    crra_utility + 1/2 * crra_utility_second_derivative * variance

  as.numeric(expected_utility)
}

get_allocations_taxable <- function(params) {
  params[1:(length(params)/2)]
}

get_allocations_taxadvantaged <- function(params) {
  params[(length(params)/2 + 1):length(params)]
}

get_taxable_indices <- function(params) {
  seq_len(length(params) / 2)
}

get_taxadvantaged_indices <- function(params) {
  seq((length(params) / 2 + 1), length(params))
}

calc_mvo_portfolio_expected_return <- function(
  params,
  expected_returns,
  ...
) {

  allocations <- params

  t(allocations) %*% expected_returns
}

calc_joint_portfolio_expected_return <- function(
  params,
  expected_returns,
  tax_matrix,
  ...
) {

  allocations_taxable       <- get_allocations_taxable(params)
  allocations_taxadvantaged <- get_allocations_taxadvantaged(params)

  t(tax_matrix %*% allocations_taxable + allocations_taxadvantaged) %*% 
  expected_returns
}

calc_joint_networth_portfolio_expected_return <- function(
  params,
  expected_returns,
  tax_matrix,
  human_capital_weights,
  liabilities_weights,
  financial_wealth,
  human_capital,
  liabilities,
  nondiscretionary_consumption,
  discretionary_consumption = 0,
  life_insurance_premium,
  income,
  ...
) {

  
  networth_fractions <- calc_networth_fractions(
    financial_wealth             = financial_wealth,
    human_capital                = human_capital,
    liabilities                  = liabilities,
    nondiscretionary_consumption = nondiscretionary_consumption,
    discretionary_consumption    = discretionary_consumption,
    income                       = income,
    life_insurance_premium       = life_insurance_premium
  )

  human_capital_frac    <- networth_fractions$human_capital
  financial_wealth_frac <- networth_fractions$financial_wealth
  liabilities_frac      <- networth_fractions$liabilities

  allocations_taxable       <- get_allocations_taxable(params)
  allocations_taxadvantaged <- get_allocations_taxadvantaged(params)

  expected_return <- 
    financial_wealth_frac * 
    t(tax_matrix %*% allocations_taxable + allocations_taxadvantaged) %*%
    expected_returns + 
    human_capital_frac * t(human_capital_weights) %*% expected_returns - 
    liabilities_frac * t(liabilities_weights) %*% expected_returns

  as.numeric(expected_return)
}

calc_mvo_portfolio_variance <- function(
  params,
  covariance_matrix,
  ...
) {

  allocations <- params
  
  t(allocations) %*% covariance_matrix %*% allocations
}

calc_joint_portfolio_variance <- function(
  params,
  covariance_matrix,
  tax_matrix,
  ...
) {
  
  allocations_taxable       <- get_allocations_taxable(params)
  allocations_taxadvantaged <- get_allocations_taxadvantaged(params)

  t(tax_matrix %*% allocations_taxable + allocations_taxadvantaged) %*% 
    covariance_matrix %*% 
    (tax_matrix %*% allocations_taxable + allocations_taxadvantaged)
}

calc_joint_networth_portfolio_variance <- function(
  params,
  covariance_matrix,
  tax_matrix,
  human_capital_weights,
  liabilities_weights,
  financial_wealth,
  human_capital,
  liabilities,
  nondiscretionary_consumption,
  discretionary_consumption = 0,
  life_insurance_premium,
  income,
  ...
) {
  
  networth_fractions <- calc_networth_fractions(
    financial_wealth             = financial_wealth,
    human_capital                = human_capital,
    liabilities                  = liabilities,
    nondiscretionary_consumption = nondiscretionary_consumption,
    discretionary_consumption    = discretionary_consumption,
    income                       = income,
    life_insurance_premium       = life_insurance_premium
  )

  human_capital_frac    <- networth_fractions$human_capital
  financial_wealth_frac <- networth_fractions$financial_wealth
  liabilities_frac      <- networth_fractions$liabilities


  allocations_taxable       <- get_allocations_taxable(params)
  allocations_taxadvantaged <- get_allocations_taxadvantaged(params)

  variance <- 
    financial_wealth_frac^2 * (
      t(tax_matrix %*% allocations_taxable + allocations_taxadvantaged) %*%
      covariance_matrix %*%
      (tax_matrix %*% allocations_taxable + allocations_taxadvantaged)
    ) + 
    2 * financial_wealth_frac * t(
      human_capital_frac * (
        covariance_matrix %*% human_capital_weights
      ) - 
      liabilities_frac * (
        covariance_matrix %*% liabilities_weights
      )
    ) %*% (
      tax_matrix %*% allocations_taxable + allocations_taxadvantaged
    ) +
    human_capital_frac^2 * (
      t(human_capital_weights) %*% covariance_matrix %*% human_capital_weights
    ) +
    liabilities_frac^2 * (
      t(liabilities_weights) %*% covariance_matrix %*% liabilities_weights
    ) - 
    2 * human_capital_frac * liabilities_frac * (
      t(human_capital_weights) %*% covariance_matrix %*% liabilities_weights
    )

  as.numeric(variance)
}

calc_networth_fractions <- function(
  financial_wealth,
  human_capital,
  liabilities,
  nondiscretionary_consumption,
  discretionary_consumption,
  income,
  life_insurance_premium
) {

  financial_wealth_prime <- 
    financial_wealth + income - discretionary_consumption - 
      nondiscretionary_consumption - life_insurance_premium
    
  human_capital_prime <- 
    human_capital - income
  
  liabilities_prime <- 
    liabilities - nondiscretionary_consumption - life_insurance_premium
    
  net_worth_prime <- 
    financial_wealth_prime + human_capital_prime - liabilities_prime

  financial_wealth_frac <- financial_wealth_prime / net_worth_prime
  human_capital_frac    <- human_capital_prime / net_worth_prime
  liabilities_frac      <- liabilities_prime / net_worth_prime

  list(
    financial_wealth = financial_wealth_frac,
    human_capital    = human_capital_frac,
    liabilities      = liabilities_frac
  )
}
