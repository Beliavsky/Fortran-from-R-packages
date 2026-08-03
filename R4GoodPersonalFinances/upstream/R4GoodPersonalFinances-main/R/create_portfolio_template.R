#' Create Portfolio Template
#' 
#' @description
#' Creates a template for default portfolio with 
#' two asset classes:
#' * `GlobalStocksIndexFund`
#' * `InflationProtectedBonds`
#' 
#' @details
#' The template is used as a starting point for creating a portfolio.
#' The asset classes have some reasonable default values of 
#' expected returns and standard deviations of returns.
#' The template assumes no correlations between asset classes in 
#' the `correlations` matrix.
#' Please check and update the template assumptions if necessary.
#' 
#' The nested `pretax` columns contain default values for
#' parameters needed for calculating effective tax rates.
#' The template assumes only capital gains tax is paid.
#' Please customise this template to your individual situation.
#' 
#' The `accounts` nested columns have zero values for all assets 
#' by default in both taxable and tax-advantaged accounts.
#' The template assumes that there is currently no financial 
#' wealth allocated to those accounts.
#' Please customise this template to your individual situation.
#' 
#' The `weights` nested columns define weights of assets in 
#' portfolios representative of the household human capital and liabilities.
#' The template assumes equal weights for all assets for both portfolios.
#' Please customise this template to your individual situation.
#' 
#' @seealso
#' Possible sources of market assumptions:
#' * \url{https://elmwealth.com/capital-market-assumptions/}  
#' * \url{https://www.obligacjeskarbowe.pl/oferta-obligacji/obligacje-10-letnie-edo/}
#' * \url{https://www.msci.com/indexes/index/664204}
#' * (PDF) https://research.ftserussell.com/Analytics/FactSheets/Home/DownloadSingleIssue?issueName=AWORLDS&isManual=False
#' 
#' @return A nested `tibble` of class 'Portfolio' with columns:
#' * `name`
#' * `expected_return`
#' * `standard_deviation`
#' * `accounts`
#'   * `taxable`
#'   * `taxadvantaged`
#' * `weights`
#'   * `human_capital`
#'   * `liabilities`
#' * `correlations`
#' * `pretax`
#'   * `turnover`
#'   * `income_qualified`
#'   * `capital_gains_long_term`
#'   * `income`
#'   * `capital_gains`
#'   * `cost_basis`
#' 
#' @examples
#'   portfolio <- create_portfolio_template()
#'   portfolio$accounts$taxable <- c(10000, 30000)
#'   portfolio
#' @export
create_portfolio_template <- function() {

    portfolio <- 
      dplyr::tribble(
        ~name,                     ~expected_return, ~standard_deviation,
        "GlobalStocksIndexFund",   0.0461,           0.15,
        "InflationProtectedBonds", 0.0200,           0,
      )
    
    correlations <- diag(1, nrow = nrow(portfolio))
    rownames(correlations) <- portfolio$name
    colnames(correlations) <- portfolio$name
  
    pretax <- 
      dplyr::tibble(
        turnover                = 0.04,
        income_qualified        = 0,
        capital_gains_long_term = 1,
        income                  = 0,
        capital_gains           = portfolio$expected_return,
        cost_basis              = 1 / ((1 + portfolio$expected_return) ^ 10)
      )
  
    portfolio <- 
      portfolio |> 
      dplyr::mutate(
        accounts = dplyr::tibble(
          taxable       = rep(0, NROW(portfolio)),
          taxadvantaged = rep(0, NROW(portfolio))
        ),
        weights = dplyr::tibble(
          human_capital = 1 / NROW(portfolio),
          liabilities   = 1 / NROW(portfolio)
        ),
        correlations = correlations,
        pretax = pretax
      )
  
    class(portfolio) <- c("Portfolio", class(portfolio))
    portfolio
}
  
#' @export
print.Portfolio <- function(x, ...) {

  name <- expected_return <- standard_deviation <- initital_value <- 
    investment_years <- preliquidation_value <- capital_gain_tax_paid <-
    postliquidation_value <- initial_value <- NULL

  args <- list(...)
  if (!is.null(args$width)) {
    class(x) <- setdiff(class(x), "Portfolio")
    print(x, ...)
    return(invisible(x))
  }

  if (is.null(args$currency)) {
    currency <- ""
  } else {
    currency <- paste0(" ", args$currency)
  }
  
  cli::cli_h1("Portfolio")

  cli::cli_h2("Market assumptions")

  market_assumptions <- 
    x |> 
    dplyr::select(name, expected_return, standard_deviation) 
  class(market_assumptions) <- setdiff(class(market_assumptions), "Portfolio")

  correlations <- x$correlations
  if (is.null(rownames(correlations))) {
    rownames(correlations) <- x$name
  }
  if (is.null(colnames(correlations))) {
    colnames(correlations) <- x$name
  }

  cli::cli_h3("Expected real returns:")
  print(market_assumptions)
  cli::cli_h3("Correlation matrix:")
  print(correlations)

  cli::cli_h2("Weights")

  weights <- as.matrix(x$weights)
  if (is.null(rownames(weights))) {
    rownames(weights) <- x$name
  }
  print(weights)

  cat("\n")
  cli::cli_h2("Accounts")

  accounts <- x$accounts
  accounts$total <- rowSums(accounts)
  accounts <- 
    rbind(
      accounts,
      colSums(accounts)
    )
  accounts <- 
    dplyr::bind_cols(
      name = c(x$name, "total"),
      accounts
    )
  print(accounts)

  cat("\n")
  cli::cli_h2("Pre-tax")
  pretax <- x$pretax
  pretax <- 
    dplyr::bind_cols(
      name = x$name,
      pretax
    )
  print(pretax, width = Inf)

  cat("\n")
  cli::cli_h2("After-tax")

  if (!"aftertax" %in% names(x)) {
    cli::cli_alert_warning(cli::col_yellow(
      "After-tax information is not available yet."
    ))
    cli::cli_alert_info(
      "Use {.code calc_effective_tax_rate()} to calculate it."
    )
  } else {
    aftertax <- x$aftertax
    aftertax <- 
      dplyr::bind_cols(
        name = x$name,
        aftertax
      ) |> 
      dplyr::select(
        -initial_value, 
        -investment_years, 
        -preliquidation_value, 
        -capital_gain_tax_paid, 
        -postliquidation_value
      )
    print(aftertax, width = Inf)
  }

  cat("\n")
  cli::cli_h2("Allocation")

  if (!"allocations" %in% names(x)) {
    cli::cli_alert_warning(cli::col_yellow(
      "Allocation information is not available yet."
    ))
    cli::cli_alert_info(
      "Use {.code calc_optimal_asset_allocation()} to calculate it."
    )
  } else {
    allocations <- x$allocations
    print(allocations, width = Inf)
  }

  invisible(x)
}
