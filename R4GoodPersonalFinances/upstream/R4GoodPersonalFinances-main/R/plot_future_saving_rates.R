#' Plotting future saving rates
#' 
#' This function plots the future saving rates from a scenario object.
#' 
#' @inheritParams plot_expected_allocation
#' @inheritParams plot_future_income
#' @param aggregation_function A function used to aggregate the saving rates
#' for multiple Monte Carlo samples. Default is `median`.
#' If `NULL`, no aggregation is performed.
#' 
#' 
#' @returns A [ggplot2::ggplot()] object. 
#' @examples 
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
#' plot_future_saving_rates(scenario)
#' @export
plot_future_saving_rates <- function(
  scenario,
  aggregation_function = stats::median,
  y_limits             = c(NA, NA)
) {

  index <- total_income <- total_spending <- savings <- saving_rate <- 
    is_mc_sample <- NULL

  colors <- PrettyCols::prettycols("Bold")

  data_to_plot <- 
    scenario |> 
    dplyr::mutate(is_mc_sample = sample > 0) 
  
  if (!is.null(aggregation_function)) {

    data_to_plot <- 
      data_to_plot |>
      dplyr::group_by(index, is_mc_sample) |>
      dplyr::summarise(
        saving_rate = aggregation_function(saving_rate)
      ) |> 
      dplyr::ungroup() 
    
    
    aggregation_function_name <- deparse(substitute(aggregation_function))
  }
  
  data_to_plot |> 
    ggplot2::ggplot(
      ggplot2::aes(
        x      = index, 
        y      = saving_rate,
        color  = saving_rate,
        shape = is_mc_sample
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(
      low     = colors[5],
      high    = colors[4],
      limits  = c(0, max(data_to_plot$saving_rate, na.rm = TRUE)),
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      limits = c(0, NA)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max(scenario$index), by = 10),
      labels = function(breaks) paste_labels(breaks, scenario = scenario)
    ) +
    ggplot2::labs(
      title = "Future Saving Rates",
      subtitle = glue::glue(paste0(
        paste_scenario_id(scenario),
        "Saving rate at year 0 is <strong>{format_percent(data_to_plot$saving_rate[1])}</strong>."
      )),
      caption = glue::glue(
        ifelse(
          max(scenario$sample) > 0,
          "Expected saving rates (round points) and from <strong>{max(scenario$sample)}</strong> Monte Carlo samples (triangles).",
          ""
        ),
        ifelse(
          max(scenario$sample) > 0 & !is.null(aggregation_function),
          paste0(
            "<br>Aggregated by function: <strong>", 
            aggregation_function_name, 
            "</strong>"
          ),
          ""
        )
      ),
      x = "Year Index",
      y = glue::glue("Saving rate"),
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position  = "none",
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    ) + 
    ggplot2::coord_cartesian(ylim = c(y_limits[1], y_limits[2])) 
}

