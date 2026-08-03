#' Calculate Effective Tax Rate
#' 
#' @inheritParams simulate_scenario
#' @param tax_rate_ltcg A numeric. Tax rate for long-term capital gains.
#' @param tax_rate_ordinary_income A numeric. Tax rate for ordinary income.
#' 
#' @return A `portfolio` object augmented with nested columns with 
#' effective tax rates calculations.
#' @examples
#'  portfolio <- create_portfolio_template()
#'  portfolio$accounts$taxable <- c(10000, 30000)
#'  portfolio <- 
#'    calc_effective_tax_rate(
#'      portfolio,
#'      tax_rate_ltcg = 0.20, 
#'      tax_rate_ordinary_income = 0.40
#'    )
#'  portfolio$aftertax$effective_tax_rate 
#' @export
calc_effective_tax_rate <- function(
  portfolio,
  tax_rate_ltcg,
  tax_rate_ordinary_income
) {
  income_qualified <- capital_gains_long_term <- blended_tax_rate_income <- blended_tax_rate_capital_gains <- income <- capital_gains <- turnover <- cost_basis <- initial_value <- preliquidation_aftertax_expected_return <- investment_years <- preliquidation_value <- capital_gain_taxed <- capital_gain_tax_paid <- postliquidation_value <- postliquidation_aftertax_expected_return <- effective_tax_rate <- NULL

  aftertax <- 
    portfolio$pretax |> 
    dplyr::transmute(
      blended_tax_rate_income = 
        income_qualified * tax_rate_ltcg + (1 - income_qualified) * tax_rate_ordinary_income,
      blended_tax_rate_capital_gains =
        capital_gains_long_term * tax_rate_ltcg + 
          (1 - capital_gains_long_term) * tax_rate_ordinary_income,
      preliquidation_aftertax_expected_return = 
        (1 - blended_tax_rate_income) * income + 
          capital_gains - 
          turnover * (1 + capital_gains - cost_basis) * blended_tax_rate_capital_gains,
      initial_value        = rep(1000, NROW(portfolio)),
      investment_years     = rep(20,   NROW(portfolio)),
      preliquidation_value = 
        initial_value * 
          (1 + preliquidation_aftertax_expected_return)^investment_years,
      capital_gain_taxed = 
        (capital_gains / portfolio$expected_return) * (1 - turnover),
      capital_gain_tax_paid = 
        (preliquidation_value - initial_value) * 
          capital_gain_taxed * tax_rate_ltcg,
      postliquidation_value = 
        preliquidation_value - capital_gain_tax_paid,
      postliquidation_aftertax_expected_return =
        (postliquidation_value / initial_value) ^ (1 / investment_years) - 1,
      effective_tax_rate = 
        1 - (postliquidation_aftertax_expected_return  / portfolio$expected_return),
      aftertax_standard_deviation = 
        (1 - effective_tax_rate) * portfolio$standard_deviation
    )
  
  portfolio <- 
    portfolio |>
    dplyr::mutate(aftertax = aftertax)

  portfolio
}
