#' Calculate optimal asset allocation
#' 
#' @inheritParams simulate_scenario
#' @returns The `portfolio` with additional nested columns:
#' * `allocations$optimal` - optimal joint net-worth portfolio allocations
#' * `allocations$current` - current allocations
#' @examplesIf interactive() 
#' older_member <- HouseholdMember$new(
#'   name       = "older",  
#'   birth_date = "1980-02-15",
#'   mode       = 80,
#'   dispersion = 10
#' )  
#' household <- Household$new()
#' household$add_member(older_member)  
#' 
#' household$expected_income <- list(
#'   "income" = c(
#'     "members$older$age <= 65 ~ 7000 * 12"
#'   )
#' )
#' household$expected_spending <- list(
#'   "spending" = c(
#'     "TRUE ~ 5000 * 12"
#'   )
#' )
#' 
#' portfolio <- create_portfolio_template() 
#' portfolio$accounts$taxable <- c(10000, 30000)
#' 
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' 
#' portfolio <- 
#'   calc_optimal_asset_allocation(
#'    household = household,
#'    portfolio = portfolio,
#'    current_date = "2020-07-15"
#'   )
#' 
#' portfolio$allocations
#' @export
calc_optimal_asset_allocation <- function(
  household,
  portfolio,
  current_date = get_current_date()
) {
  asset_class <- NULL

  scenario <- 
    simulate_scenario(
      household    = household,
      portfolio    = portfolio,
      current_date = current_date
    )

  financial_wealth <- scenario$financial_wealth[1]
  human_capital    <- scenario$human_capital[1]
  liabilities      <- scenario$liabilities[1]
  net_worth        <- financial_wealth + human_capital - liabilities
  
  nondiscretionary_consumption <- scenario$nondiscretionary_spending[1]

  fraction_in_taxable_accounts <- 
    sum(portfolio$accounts$taxable) / sum(portfolio$accounts)
  
  discretionary_spending <- scenario$discretionary_spending[1]
  income                 <- scenario$total_income[1]
  
  optimal_joint_networth_portfolio <- calc_optimal_portfolio(
    risk_tolerance               = household$risk_tolerance,
    expected_returns             = portfolio$expected_return,
    standard_deviations          = portfolio$standard_deviation,
    correlations                 = portfolio$correlations,
    effective_tax_rates          = portfolio$aftertax$effective_tax_rate,
    in_taxable_accounts          = fraction_in_taxable_accounts,
    financial_wealth             = financial_wealth,
    human_capital                = human_capital,
    liabilities                  = liabilities,
    nondiscretionary_consumption = nondiscretionary_consumption,
    discretionary_consumption    = discretionary_spending,
    income                       = income,
    life_insurance_premium       = 0,
    human_capital_weights        = portfolio$weights$human_capital,
    liabilities_weights          = portfolio$weights$liabilities,
    asset_names                  = portfolio$name
  ) 

  current_allocation <- 
    portfolio$accounts / sum(portfolio$accounts)

  current_allocation <- 
    current_allocation |> 
    dplyr::mutate(total = rowSums(current_allocation))

  portfolio <-
    portfolio |>
    dplyr::mutate(
      allocations = dplyr::tibble(
        name = optimal_joint_networth_portfolio$allocations$asset_class,
        optimal = 
          optimal_joint_networth_portfolio$allocations |> 
            dplyr::select(-asset_class),
        current = 
          current_allocation |> 
          dplyr::as_tibble()
      )
    )
  
  portfolio
}
