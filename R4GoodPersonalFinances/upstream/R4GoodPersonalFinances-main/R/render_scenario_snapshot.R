#' Rendering a scenario snapshot
#' 
#' @inheritParams plot_expected_allocation
#' @param index The index of the scenario year to render.
#' By default, it is 0, which corresponds to the current year.
#' @param currency The currency symbol to use as a suffix.
#' @param big_mark The character to use as a big mark.
#' It separates thousands.
#' @returns A [gt::gt()] object.
#' @examples
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
#'     "members$older$age <= 65 ~ 9000 * 12"
#'   )
#' )
#' household$expected_spending <- list(
#'   "spending" = c(
#'     "members$older$age <= 65 ~ 5000 * 12",
#'     "TRUE ~ 4000 * 12"
#'   )
#' )
#' 
#' portfolio <- create_portfolio_template() 
#' portfolio$accounts$taxable <- c(10000, 30000)
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' 
#' scenario <- 
#'   simulate_scenario(
#'    household = household,
#'    portfolio = portfolio,
#'    current_date = "2020-07-15"
#'   )
#' render_scenario_snapshot(scenario)
#' @export
render_scenario_snapshot <- function(
  scenario,
  index    = 0, 
  currency = "", 
  big_mark = " "
) {

  value <- NULL

  scenario_at_index <- 
    scenario |> 
    dplyr::filter(index == !!index)

  currency <- paste0(" ", currency)

  tibble::tribble(
    ~name, ~value,
    "Financial wealth", 
      format_currency(
        scenario_at_index$financial_wealth, 
        suffix = currency, 
        big.mark = big_mark
      ),
    "Human capital", 
      format_currency(
        scenario_at_index$human_capital, 
        suffix = currency, 
        big.mark = big_mark
      ),
    "Liabilities", 
      format_currency(
        scenario_at_index$liabilities, 
        suffix = currency, 
        big.mark = big_mark
      ),
    "Net Worth", 
      format_currency(
        scenario_at_index$net_worth, 
        suffix = currency, 
        big.mark = big_mark
      ),
    "Income", 
    format_currency(
      scenario_at_index$total_income / 12, 
      accuracy = 1, 
      suffix = currency, big.mark = big_mark
    ),
    "Spending", 
    format_currency(
      scenario_at_index$total_spending / 12, 
      accuracy = 1, 
      suffix = currency, big.mark = big_mark
    ),
    "Non-discretionary spending", 
      format_currency(
        scenario_at_index$nondiscretionary_spending / 12, 
        accuracy = 1, 
        suffix = currency, big.mark = big_mark
      ),
    "Discretionary spending", 
      format_currency(
        scenario_at_index$discretionary_spending / 12, 
        accuracy = 1, 
        suffix = currency, big.mark = big_mark
      ),
    "Savings", format_currency(
      scenario_at_index$savings / 12, 
      accuracy = 1, 
      suffix = currency, big.mark = big_mark
    ),
    "Saving rate", format_percent(scenario_at_index$saving_rate),
  ) |> 
    tibble::column_to_rownames(var = "name") |> 
    gt::gt(rownames_to_stub = TRUE) |> 
    gt::cols_label(value = "") |>
    gt::tab_header(
      subtitle = gt::md(glue::glue("Scenario: **{scenario_at_index$scenario_id}**")),
      title = gt::md(glue::glue(
        "Scenario Summary ({scenario_at_index$year})"
      ))
    ) |>
    gt::tab_source_note(
      source_note = gt::md(glue::glue("*Date: {scenario_at_index$date}*"))
    ) |> 
    gt::tab_row_group(
      label = gt::md("*Balance sheet*"),
      rows = c("Financial wealth", "Human capital", "Liabilities", "Net Worth")
    ) |> 
    gt::tab_row_group(
      label = gt::md("*Expected cashflow (monthly)*"),
      rows = c(
        "Income", "Spending", 
        c("Savings", "Saving rate"),
        c("Non-discretionary spending", "Discretionary spending")
      ) 
    ) |> 
    gt::data_color(
      rows = dplyr::starts_with(c("Discretionary")),
      palette = PrettyCols::prettycols("Bold")[c(1)]
    ) |> 
    gt::data_color(
      rows = dplyr::starts_with(c("Savings")),
      palette = PrettyCols::prettycols("Bold")[c(4)]
    ) |> 
    gt::cols_align(columns = value, align = "right")
}
