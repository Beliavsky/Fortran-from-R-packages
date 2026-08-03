#' Plot scenarios metrics
#' 
#' @description
#' The plot allows to compare metrics for multiple scenarios.
#' 
#' If scenarios are simulated without Monte Carlo samples,
#' so they are based only on expected returns of portfolio,
#' two metrics are available for each scenario:
#' * constant discretionary spending - certainty equivalent constant
#' level of consumption that would result in the same lifetime utility 
#' as a given series of future consumption in a given scenario
#' (the higher, the better).
#' * utility of discretionary spending - normalized
#' to minimum and maximum values of constant discretionary spending
#' (the higher, the better).
#' 
#' If scenarios are simulated with additional Monte Carlo samples,
#' there are four more metrics available per scenario:
#' * constant discretionary spending (for Monte Carlo samples),
#' * normalized median utility of discretionary spending 
#' (for Monte Carlo samples),
#' * median of missing funds that need additional income 
#' or additional savings at the expense of non-discretionary spending,
#' (of yearly averages of Monte Carlo samples),
#' * median of discretionary spending 
#' (of yearly averages of Monte Carlo samples).
#' 
#' 
#' @param scenarios A `tibble` with nested columns - 
#' the result of [simulate_scenarios()].
#' @param period A character. The amounts can be shown
#' as yearly values (default) or averaged per month values.
#' @returns A [ggplot2::ggplot()] object. 
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
#'     "is_not_on('older', 'retirement') ~ 7000 * 12"
#'   )
#' )
#' household$expected_spending <- list(
#'   "spending" = c(
#'     "TRUE ~ 4000 * 12"
#'   )
#' )
#' 
#' portfolio <- create_portfolio_template() 
#' portfolio$accounts$taxable <- c(100000, 300000)
#' portfolio <- 
#'   portfolio |> 
#'   calc_effective_tax_rate(
#'     tax_rate_ltcg = 0.20, 
#'     tax_rate_ordinary_income = 0.40
#'   )
#' 
#' start_ages <- c(60, 65, 75)
#' scenarios_parameters <- 
#'   tibble::tibble(
#'     member    = "older",
#'     event      = "retirement",
#'     start_age = start_ages,
#'     years     = Inf,
#'     end_age   = Inf
#'    ) |> 
#'   dplyr::mutate(scenario_id = start_age) |> 
#'   tidyr::nest(events = -scenario_id)
#' 
#' scenarios <- 
#'   simulate_scenarios(
#'     scenarios_parameters = scenarios_parameters,
#'     household            = household,
#'     portfolio            = portfolio,
#'     maxeval              = 100,
#'     current_date         = "2020-07-15"
#'   )
#' 
#' plot_scenarios(scenarios, "monthly")
#' @export
plot_scenarios <- function(
  scenarios,
  period = c("yearly", "monthly")
) {

  scenario_id <- discretionary_spending_utility_weighted <- survival_prob <- 
    time_value_discount <- smooth_consumption_preference <- utility_expected <- 
    constant_expected <- utility <- constant <- 
    positive_discretionary_spending <- negative_discretionary_spending <- 
    risk_tolerance <- discretionary_spending <- utility_normalized_expected <- 
    utility_normalized <- metric <- value <- NULL

  period        <- rlang::arg_match(period)
  period_factor <- if (period == "yearly") 1 else 12
    
  expected_returns_scenario <- 
    scenarios |> 
    dplyr::filter(sample == 0) |> 
    dplyr::group_by(scenario_id) |> 
    dplyr::summarise(
      utility_expected  = sum(discretionary_spending_utility_weighted),
      constant_expected = 
        calc_inverse_utility(
          utility = 
            sum(discretionary_spending_utility_weighted) / 
            sum(survival_prob * time_value_discount),
          parameter = unique(smooth_consumption_preference)
        ) / period_factor
    )
  
  monte_carlo_scenarios <-
    scenarios |> 
    dplyr::filter(sample != 0) 
  
  if (NROW(monte_carlo_scenarios) > 0) {
  
    monte_carlo_scenarios <-
      monte_carlo_scenarios |> 
      dplyr::mutate(
        negative_discretionary_spending = 
          dplyr::if_else(
            discretionary_spending < 0,
            abs(discretionary_spending),
            0
          ),
        positive_discretionary_spending = 
          dplyr::if_else(
            discretionary_spending >= 0,
            discretionary_spending,
            0
          )
        ) |> 
      dplyr::group_by(scenario_id, sample) |> 
      dplyr::summarise(
        utility = 
          sum(discretionary_spending_utility_weighted),
        constant = 
          calc_inverse_utility(
            utility = 
              sum(discretionary_spending_utility_weighted) / 
              sum(survival_prob * time_value_discount),
            parameter = unique(smooth_consumption_preference)
          ),
        positive_discretionary_spending = mean(positive_discretionary_spending),
        negative_discretionary_spending = mean(negative_discretionary_spending),
        risk_tolerance = unique(risk_tolerance)
      ) |> 
      dplyr::group_by(scenario_id) |> 
      dplyr::summarise(
        utility = stats::median(utility),
        constant = 
          calc_inverse_utility(
            mean(
              calc_utility(
              constant,
              parameter = unique(risk_tolerance)
              )
            ),
            parameter = unique(risk_tolerance)
          ) / period_factor,
        positive_discretionary_spending = 
          stats::median(positive_discretionary_spending) / period_factor,
        negative_discretionary_spending = 
          stats::median(negative_discretionary_spending) / period_factor
      ) 
  }
  
  expected_returns_scenario <- 
    expected_returns_scenario |> 
      dplyr::mutate(
        utility_normalized_expected = 
          normalize(
            utility_expected, 
            min_val = min(
              constant_expected, 
              ifelse(
                NROW(monte_carlo_scenarios) == 0,
                NA, 
                monte_carlo_scenarios$constant
              ),
              na.rm = TRUE
            ), 
            max_val = max(
              constant_expected, 
              ifelse(
                NROW(monte_carlo_scenarios) == 0,
                NA, 
                monte_carlo_scenarios$constant
              ),
              na.rm = TRUE
            )
          )
      )
  
  if (NROW(monte_carlo_scenarios) > 0) {
    monte_carlo_scenarios <- 
      monte_carlo_scenarios |> 
      dplyr::mutate(
        utility_normalized = 
          normalize(
            utility, 
            min_val = min(
              constant, 
              expected_returns_scenario$constant_expected,
              na.rm = TRUE
            ), 
            max_val = max(
              constant, 
              expected_returns_scenario$constant_expected,
              na.rm = TRUE
            )
          )
      )
  }

  ordered_scenario_levels <- 
    expected_returns_scenario |>
    dplyr::pull(scenario_id) |> 
    unique()

  expected_returns_scenario_long <- 
    expected_returns_scenario |> 
    dplyr::select(-utility_expected) |> 
    tidyr::pivot_longer(
      cols      = -"scenario_id",
      names_to  = "metric",
      values_to = "value"
    ) |> 
      dplyr::mutate(
        scenario_id = factor(scenario_id, levels = ordered_scenario_levels)
      ) 
    
  if (NROW(monte_carlo_scenarios) > 0) {
    
    monte_carlo_scenarios_long <- 
      monte_carlo_scenarios |> 
      dplyr::select(-utility) |> 
      tidyr::pivot_longer(
        cols      = -"scenario_id",
        names_to  = "metric",
        values_to = "value"
      ) |> 
      dplyr::mutate(
        scenario_id = factor(scenario_id, levels = ordered_scenario_levels)
      ) 
  }
    
  colors <- 
    grDevices::colorRampPalette(
      rev(PrettyCols::prettycols("Bold"))
    )(4 + 2)

  to_plot <- 
    expected_returns_scenario_long |>
    ggplot2::ggplot(
      ggplot2::aes(
        x     = scenario_id, 
        y     = value,
        color = metric,
        group = metric
      )
    ) +
    ggplot2::geom_line(
      data     = expected_returns_scenario_long,
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      data = expected_returns_scenario_long,
      size = 2
    ) +
    ggplot2::geom_label(
      data = 
        expected_returns_scenario_long |> 
        dplyr::filter(metric == "constant_expected"),
      ggplot2::aes(
        label = 
          ifelse(
            round(value / 1000, 1)  == 0, 
            "0k",
            paste0(round(value / 1000, 1), "k")
          )
      ),
      nudge_x     = -0.2,
      na.rm       = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = format_currency,
      breaks = scales::breaks_extended(n = 10)
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      labels = c(
        "negative_discretionary_spending" = 
          "median of means\nof missing funds\nin Monte Carlo samples",
        "positive_discretionary_spending" = 
          "median of means\nof discretionary spending\nin Monte Carlo samples",
        "constant" = 
          "constant (certainty equivalent)\ndiscretionary spending\n in Monte Carlo samples",
        "constant_expected" = 
          "constant (certainty equivalent)\ndiscretionary spending\nbased on expected returns",
        "utility_normalized_expected" = 
          "normalized utility\nof discretionary spending\nbased on expected returns",
        "utility_normalized" = 
          "normalized median utility\nof discretionary spending\nin Monte Carlo samples"
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title    = glue::glue("Comparison of Scenario Metrics"),
      subtitle = glue::glue("Spending period: <strong>{period}</strong>"),
      x        = "Scenario",
      y        = glue::glue("Amount ({period})"),
    ) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(
        angle = 70, 
        hjust = 1
      ),
      plot.caption = 
        ggtext::element_markdown(
          color = "grey60", 
          size  = 10
        ),
      plot.subtitle = ggtext::element_markdown(color = "grey60")
    )

    if (NROW(monte_carlo_scenarios) > 0) {
      
      to_plot <- 
        to_plot +
        ggplot2::geom_line(
          data    = monte_carlo_scenarios_long,
          linetype = "dotted"
        ) +
        ggplot2::geom_point(
          data = monte_carlo_scenarios_long,
          size = 2
        )
    }
    
  to_plot
}
