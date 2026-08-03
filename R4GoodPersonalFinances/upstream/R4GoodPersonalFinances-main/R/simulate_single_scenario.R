simulate_single_scenario <- function(
  household,
  portfolio,
  scenario_id    = "default",
  current_date   = get_current_date(),
  random_returns = FALSE,
  seed           = NULL,
  debug          = FALSE,
  ...
) {

  income <- spending <- total_income <- nondiscretionary_spending <- 
    consumption_impatience_preference <- index <- 
    smooth_consumption_preference <- survival_prob <- time_value_discount <- 
    discretionary_spending_utility <- total_spending <- 
    savings <- saving_rate <- NULL

  current_date <- lubridate::as_date(current_date)

  timeline <- 
    generate_household_timeline(
      household    = household, 
      current_date = current_date
    ) 
  
  income_streams <- 
    generate_cashflow_streams(
      timeline = timeline,  
      triggers = household$expected_income
    )
  spending_streams <- 
    generate_cashflow_streams(
      timeline = timeline,  
      triggers = household$expected_spending
    )
  
  scenario <- 
    timeline |> 
    dplyr::mutate(
      "income"                  = income_streams,
      total_income              = rowSums(income),
      "spending"                = spending_streams,
      nondiscretionary_spending = rowSums(spending)
    ) 

  fraction_in_taxable_accounts <- 
    sum(portfolio$accounts$taxable) / sum(portfolio$accounts)
  
  financial_wealth <- sum(portfolio$accounts)
  
  weights <- list()

  if (financial_wealth > 0) {

    weights$taxable <- 
      portfolio$accounts$taxable / financial_wealth
    weights$taxadvantaged <- 
      portfolio$accounts$taxadvantaged / financial_wealth
    weights$financial_wealth <- 
      weights$taxable + weights$taxadvantaged

  } else {

    weights$taxable          <- rep(0, NROW(portfolio))
    weights$taxadvantaged    <- rep(0, NROW(portfolio))
    weights$financial_wealth <- rep(0, NROW(portfolio))
  }

  human_capital_discount_rate <- 
    calc_portfolio_expected_return(
      weights = portfolio$weights$human_capital,
      returns = portfolio$expected_return
    )
  
  liabilities_discount_rate <- 
    calc_portfolio_expected_return(
      weights = portfolio$weights$liabilities,
      returns = portfolio$expected_return
    )

  if (random_returns) {

      returns <- 
        generate_random_returns(
          portfolio = portfolio, 
          n         = NROW(scenario),
          seed      = seed
        )
    
  } else {

    returns <- t(portfolio$expected_return) 
    colnames(returns) <-portfolio$name
    returns <- dplyr::as_tibble(returns) 
    returns <- returns[rep(1, NROW(scenario)), ]
  }
  
  scenario <- 
    scenario |>
    dplyr::mutate(
      human_capital = calc_present_value(
        cashflow      = total_income,
        discount_rate = human_capital_discount_rate
      ),
      liabilities = calc_present_value(
        cashflow      = nondiscretionary_spending,
        discount_rate = liabilities_discount_rate
      ),
      portfolio = dplyr::tibble(
        allocation = rep(list(NA_real_), NROW(scenario)),
        returns    = returns,
      ),
      financial_wealth       = NA_real_,
      net_worth              = NA_real_,
      discretionary_spending = NA_real_,
      total_spending         = NA_real_,
      financial_wealth_end   = NA_real_,
      risk_tolerance         = household$risk_tolerance,
      smooth_consumption_preference = 
        household$smooth_consumption_preference,
      consumption_impatience_preference = 
        household$consumption_impatience_preference
    ) 
  
  household_survival <- household$calc_survival(current_date = current_date)
  household_min_age  <- household$get_min_age(  current_date = current_date)

  max_age             <- household_min_age + max(scenario$index)
  gompertz_mode       <- household_survival$mode
  gompertz_dispersion <- household_survival$dispersion

  n_rows <- NROW(scenario)

  for (i in seq_len(n_rows)) {

    if (debug) {
      cli::cli_alert(
        "Index year: {i}/{n_rows}." 
      )
    }

    if (i == 1) {

      scenario[i, ]$financial_wealth <- financial_wealth

      scenario[i, ]$portfolio$allocation <- list(
        dplyr::tibble(
          asset         = portfolio$name,
          taxable       = weights$taxable,
          taxadvantaged = weights$taxadvantaged,
          total         = weights$financial_wealth
        )
      )
      initial_allocation <- c(
        scenario[i, ]$portfolio$allocation[[1]]$taxable,
        scenario[i, ]$portfolio$allocation[[1]]$taxadvantaged
      )
      
    } else {

      initial_allocation <- scenario[i - 1, ]$portfolio$allocation[[1]]
      initial_allocation <- c(
        initial_allocation$taxable,
        initial_allocation$taxadvantaged
      )
    }

    scenario[i, ]$net_worth <- 
      scenario[i, ]$financial_wealth +
      scenario[i, ]$human_capital -
      scenario[i, ]$liabilities

    discretionary_spending <- 
      calc_discretionary_spending(
        allocations_taxable          = weights$taxable,
        allocations_taxadvantaged    = weights$taxadvantaged,
        human_capital_weights        = portfolio$weights$human_capital,
        liabilities_weights          = portfolio$weights$liabilities,
        expected_returns             = portfolio$expected_return,
        standard_deviations          = portfolio$standard_deviation,
        effective_tax_rates          = portfolio$aftertax$effective_tax_rate,
        correlations                 = portfolio$correlations,
        financial_wealth             = scenario[i, ]$financial_wealth,
        human_capital                = scenario[i, ]$human_capital,
        liabilities                  = scenario[i, ]$liabilities,
        nondiscretionary_consumption = scenario[i, ]$nondiscretionary_spending,
        income                       = scenario[i, ]$total_income,
        risk_tolerance               = scenario[i, ]$risk_tolerance,
        consumption_impatience_preference = 
          scenario[i, ]$consumption_impatience_preference,
        smooth_consumption_preference = 
          scenario[i, ]$smooth_consumption_preference,
        current_age = household_min_age + scenario[i, ]$index,
        max_age                = max_age,
        gompertz_mode          = gompertz_mode,
        gompertz_dispersion    = gompertz_dispersion,
        life_insurance_premium = 0
      )
    
    if (is.infinite(discretionary_spending)) {

      if (debug) {
        cli::cli_alert_warning(
          cli::col_yellow(paste(
            "Discretionary spending is infinite.",
            "Setting it to this period net-worth."
          ))
        )
      }

      discretionary_spending <- 
        scenario[i, ]$net_worth
    }

    scenario[i, ]$discretionary_spending <- discretionary_spending

    scenario[i, ]$total_spending <- 
      scenario[i, ]$discretionary_spending +
      scenario[i, ]$nondiscretionary_spending

    financial_wealth_end <- 
      scenario[i, ]$financial_wealth +
      scenario[i, ]$total_income -
      scenario[i, ]$total_spending
      
    optimal_joint_networth_portfolio <- tryCatch({

      calc_optimal_portfolio(
        risk_tolerance               = household$risk_tolerance,
        expected_returns             = portfolio$expected_return,
        standard_deviations          = portfolio$standard_deviation,
        correlations                 = portfolio$correlations,
        effective_tax_rates          = portfolio$aftertax$effective_tax_rate,
        in_taxable_accounts          = fraction_in_taxable_accounts,
        financial_wealth             = scenario[i, ]$financial_wealth,
        human_capital                = scenario[i, ]$human_capital,
        liabilities                  = scenario[i, ]$liabilities,
        nondiscretionary_consumption = scenario[i, ]$nondiscretionary_spending,
        discretionary_consumption    = scenario[i, ]$discretionary_spending,
        income                       = scenario[i, ]$total_income,
        life_insurance_premium       = 0,
        human_capital_weights        = portfolio$weights$human_capital,
        liabilities_weights          = portfolio$weights$liabilities,
        asset_names                  = portfolio$name,
        initial_allocation           = initial_allocation,
        ...
        )
    },

        error = function(e) {
          if (debug) {
            cli::cli_alert_warning(
              cli::col_yellow(
                "Optimal allocation not found." 
              )
            )
          }
          return(NULL)
        }
    )
      
    if (
      !is.null(optimal_joint_networth_portfolio$allocations$total) && 
        i < n_rows
    ) {

      scenario[i, ]$portfolio$allocation <- list(
        optimal_joint_networth_portfolio$allocations |> 
          dplyr::rename("asset" = "asset_class")
      )

    } else if (i > 1) {

      if (debug) {
        cli::cli_alert_warning(
          cli::col_yellow(
            "Using optimal allocation from previous period ({i - 1})."
          )
        )
      }

      scenario[i, ]$portfolio$allocation <- 
        scenario[i - 1, ]$portfolio$allocation

    } else (

      return(NULL)
    )

    financial_wealth_end <- 
      sum(
        rep(financial_wealth_end, NROW(portfolio)) *
          scenario[i, ]$portfolio$allocation[[1]]$total * 
          (1 + t(scenario[i, ]$portfolio$returns))
      )
      
    scenario[i, ]$financial_wealth_end <- financial_wealth_end

    if (i < n_rows) {
      scenario[i + 1, ]$financial_wealth <- financial_wealth_end
    }
  }

  scenario |> 
    dplyr::mutate(
      time_value_discount = 
        1 / (
          (1 + consumption_impatience_preference)^(index - 0)
        ),
      discretionary_spending_utility = 
        calc_utility(
          x         = discretionary_spending, 
          parameter = smooth_consumption_preference
        ),
      discretionary_spending_utility_weighted = 
        survival_prob * time_value_discount * discretionary_spending_utility,
      scenario_id = scenario_id,
      savings     = dplyr::if_else(
        total_income > 0 & total_spending > 0,
        total_income - total_spending,
        0
      ),
      savings = dplyr::if_else(savings < 0, 0, savings),
      saving_rate = dplyr::if_else(
        total_income > 0,
        savings / total_income,
        0
      ),
      saving_rate = dplyr::if_else(saving_rate < 0, 0, saving_rate)
    ) |> 
    dplyr::select(
      scenario_id, 
      index, 
      dplyr::everything()
    )
}
