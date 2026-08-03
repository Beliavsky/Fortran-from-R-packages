#' Plot expected capital over household life cycle
#' 
#' Plots financial capital, human capital, total capital, and liabilities.
#' @inheritParams plot_expected_allocation
#' @returns A [ggplot2::ggplot()] object.
#' @examplesIf interactive()
#' older_member <- HouseholdMember$new(
#'   name       = "older",  
#'   birth_date = "2000-02-15",
#'   mode       = 80,
#'   dispersion = 10
#' )  
#' household <- Household$new()
#' household$add_member(older_member)  
#' 
#' household$expected_income <- list(
#'   "income" = c(
#'     "members$older$age <= 65 ~ 10000 * 12"
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
#' portfolio$weights$human_capital <- c(0.2, 0.8)
#' portfolio$weights$liabilities <- c(0.1, 0.9)
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg            = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' 
#' scenario <- 
#'   simulate_scenario(
#'    household    = household,
#'    portfolio    = portfolio,
#'    current_date = "2020-07-15"
#'   )
#' 
#' plot_expected_capital(scenario)
#' @export
plot_expected_capital <- function(
  scenario
) {
  index <- financial_wealth <- human_capital <- liabilities <- amount <- 
    type <- NULL

  colors <- PrettyCols::prettycols("Bold")
  color_values <- c(
    "financial_wealth" = colors[1],
    "human_capital"    = colors[4],
    "liabilities"      = colors[5],
    "total_capital"    = "gray85"
  )
  color_labels <- c(
    "financial_wealth" = "Financial capital",
    "human_capital"    = "Human capital",
    "liabilities"      = "Liabilities",
    "total_capital"    = "Total capital"
  )

  scenario |> 
    dplyr::filter(sample == 0) |>
    dplyr::select(
      index,
      financial_wealth, 
      human_capital, 
      liabilities
    ) |>
    dplyr::add_row(
      index            = max(!!scenario$index) + 1,
      financial_wealth = 0,
      human_capital    = 0,
      liabilities      = 0
    ) |> 
    dplyr::mutate(
      total_capital = financial_wealth + human_capital
    ) |> 
    tidyr::pivot_longer(
      cols      = -c("index"),
      names_to  = "type",
      values_to = "amount"
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x     = index, 
        y     = amount, 
        color = type,
        fill  = type
      )
    ) +
    ggplot2::geom_area(
      alpha    = 0.20, 
      position = "identity"
    ) +
    ggplot2::scale_y_continuous(
      labels = format_currency
    ) +
    ggplot2::scale_color_manual(
      values = color_values,
      labels = color_labels
    ) +
    ggplot2::scale_fill_manual(
      values = color_values,
      labels = color_labels
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(scenario$index), by = 10),
      labels = 
        function(breaks) paste_labels(breaks, scenario = scenario)
    ) +
    ggplot2::labs(
      title = glue::glue(
        "Expected Financial and Human Capital over Household Life Cycle"
      ),
      subtitle = paste_scenario_id(scenario),
      x = paste_year_index_axis_label(),
      y = "Amount",
    ) +
    ggplot2::theme(
      legend.position  = "right",
      legend.title     = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    )
    
}
